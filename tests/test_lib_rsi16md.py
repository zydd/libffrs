
#  test_lib_rsi16md.py
#
#  Copyright 2026 Gabriel Machado
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import pytest
import random

import ffrs
import ffrs.reference as ref
import ffrs.reference.rs as ref_rs
import ffrs.reference.ntt
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted, randbytes


ref_gf = ref.GF(65537, 1, 3)

ref_ntt = lambda w, buf: to_int_list(ffrs.reference.ntt.ntt(ref_gf, ref_gf(w), rbo_sorted(buf)))
ref_intt = lambda w, buf: rbo_sorted(to_int_list(ffrs.reference.ntt.intt(ref_gf, ref_gf(w), buf)))


class BaseTestRS:
    def test_encode(self, rs):
        buf = randbytes(rs.message_size)

        res = rs.encode(buf)

        w = ref_gf(rs.root)
        buf_ntt = ref_ntt(w, to_int_list(buf + bytearray(rs.ntt_size - len(buf)), 2))

        ecc_mix = rs._mix_ecc(buf_ntt[:rs.ecc_len])

        assert to_bytearray(ecc_mix, 2) == res

    @pytest.mark.parametrize("grace", [0, 1, 2, 3, 4])
    def test_repair_unknown_locations(self, rs, grace):
        orig = [random.randrange(0, 2**16) for _ in range(rs.block_len - rs.ecc_len)]

        block_enc = orig + to_int_list(rs.encode(to_bytearray(orig, 2)), 2)

        block_enc_err = list(block_enc)

        error_positions = dict()

        # # Always add one error to the codeword
        # if len(error_positions) < max(rs.ecc_len//2, 1):
        #     i = random.randrange(rs.block_len - rs.ecc_len, rs.block_len)
        #     error_positions[i] = random.randrange(2**16)
        #     block_enc_err[i] = rs.gf.add(block_enc_err[i], error_positions[i])

        while len(error_positions) < max(rs.ecc_len//2 - grace, 1):
            i = random.randrange(rs.block_len)
            if i in error_positions:
                continue
            error_positions[i] = random.randrange(2**16)
            block_enc_err[i] = rs.gf.add(block_enc_err[i], error_positions[i])

        # print("error_positions:", error_positions)

        block_enc_err = to_bytearray(block_enc_err, 2)
        msg_err = block_enc_err[:rs.message_size]
        ecc_err = block_enc_err[rs.message_size:]

        assert block_enc != to_int_list(msg_err + ecc_err, 2)

        # synd = to_int_list(rs.synd(msg_err, ecc_err), 2)
        # loc, ev = ref_rs.sugiyama(ref_gf, synd)
        # rec_pos, rec_val = ref_rs.forney(rs.block_len, ref_gf(rs.root), loc, ev)

        rs.repair(msg_err, ecc_err)

        assert block_enc == to_int_list(msg_err + ecc_err, 2)

    def test_repair(self, rs):
        orig = [random.randrange(0, 2**16) for _ in range(rs.block_len - rs.ecc_len)]

        block_enc = orig + to_int_list(rs.encode(to_bytearray(orig, 2)), 2)

        block_enc_err = list(block_enc)

        error_positions = dict()

        # Always add one error to the codeword
        if len(error_positions) < max(rs.ecc_len, 1):
            i = random.randrange(rs.block_len - rs.ecc_len, rs.block_len)
            error_positions[i] = random.randrange(2**16)
            block_enc_err[i] = rs.gf.add(block_enc_err[i], error_positions[i])

        while len(error_positions) < rs.ecc_len:
            i = random.randrange(rs.block_len)
            if i in error_positions:
                continue
            error_positions[i] = random.randrange(2**16)
            block_enc_err[i] = rs.gf.add(block_enc_err[i], error_positions[i])

        w = ref_gf(rs.root)
        print()
        print("errs: ", error_positions)

        block_enc_err = to_bytearray(block_enc_err, 2)
        msg_err = block_enc_err[:rs.message_size]
        ecc_err = block_enc_err[rs.message_size:]

        assert block_enc != to_int_list(msg_err + ecc_err, 2)

        rs.repair(msg_err, ecc_err, list(error_positions.keys()))

        assert block_enc == to_int_list(msg_err + ecc_err, 2)

    def test_repair_no_errors(self, rs):
        orig = [random.randrange(0, 2**16) for _ in range(rs.block_len - rs.ecc_len)]
        msg = to_bytearray(orig, 2)
        ecc = rs.encode(msg)
        msg_err = bytearray(msg)
        ecc_err = bytearray(ecc)

        rs.repair(msg_err, ecc_err, [])

        assert msg_err == msg
        assert ecc_err == ecc

    def test_repair_unknown_no_errors(self, rs):
        orig = [random.randrange(0, 2**16) for _ in range(rs.block_len - rs.ecc_len)]
        msg = to_bytearray(orig, 2)
        ecc = rs.encode(msg)
        msg_err = bytearray(msg)
        ecc_err = bytearray(ecc)

        rs.repair(msg_err, ecc_err)

        assert msg_err == msg
        assert ecc_err == ecc

    def test_encode_blocks_empty(self, rs):
        assert bytearray() == rs.encode(bytearray())

    def test_encode_blocks_single(self, rs):
        buf = randbytes(rs.message_size)

        buf_enc = rs.encode(buf)
        buf_enc_blk = rs.encode(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_size
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize("count", [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 100])
    def test_encode_blocks_multiple(self, rs, count):
        buf = randbytes(rs.message_size * count)

        buf_enc = [rs.encode(buf[i * rs.message_size:(i + 1) * rs.message_size])
                        for i in range(count)]
        buf_enc = b"".join(buf_enc)
        buf_enc_blk = rs.encode(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_size * count
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize("extra", [2, 4, 6, 8, 10, 12, 14, 128, 1500])
    def test_encode_blocks_remainder(self, rs, extra):
        count = 16 + extra // rs.message_size
        extra = (extra % rs.message_size) or 2
        buf = randbytes(rs.message_size * count + extra)

        buf_enc = [rs.encode(buf[i * rs.message_size:(i + 1) * rs.message_size])
                        for i in range(count)]

        buf_enc += [rs.encode(buf[-extra:] + bytearray(rs.message_size - extra))]
        buf_enc = b"".join(buf_enc)
        buf_enc_blk = rs.encode(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_size * (count + 1)
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize("interleave", list(range(32)) + [32, 48, 100, 256])
    def test_encode_interleaved(self, rs, interleave):
        data = list(range(interleave * rs.message_len))
        buf = to_bytearray(data, 2)

        assert rs.interleave == 1
        rs.interleave = interleave
        try:
            interleaved_chunk = rs.encode(buf)
        finally:
            rs.interleave = 1

        rs.interleave = 1
        for i in range(interleave):
            buf_enc = rs.encode(to_bytearray(data[i::interleave], 2))

            # Sequential ecc
            # assert interleaved_chunk[i * rs.ecc_size:(i + 1) * rs.ecc_size] == buf_enc

            # Interleaved ecc
            assert to_bytearray(to_int_list(interleaved_chunk, 2)[i::interleave], 2) == buf_enc


@pytest.mark.parametrize("rs", [
    ffrs.RSi16md(4, ecc_len=2),
    ffrs.RSi16md(16, ecc_len=4),
    ffrs.RSi16md(16, ecc_len=8),
    ffrs.RSi16md(128, ecc_len=2),
    ffrs.RSi16md(256, ecc_len=128),
    # ffrs.RSi16md(1024, ecc_len=128),
    # ffrs.RSi16md(4096, ecc_len=512),

    ffrs.RSi16md(4, ecc_len=2, simd_x16=False),
    ffrs.RSi16md(16, ecc_len=4, simd_x16=False),
    ffrs.RSi16md(128, ecc_len=2, simd_x16=False),
    ffrs.RSi16md(256, ecc_len=128, simd_x16=False),
    # ffrs.RSi16md(1024, ecc_len=128, simd_x16=False),
    # ffrs.RSi16md(4096, ecc_len=512, simd_x16=False),

    ffrs.RSi16md(4, ecc_len=2, simd_x16=False, simd_x8=False),
    ffrs.RSi16md(16, ecc_len=4, simd_x16=False, simd_x8=False),
    ffrs.RSi16md(128, ecc_len=2, simd_x16=False, simd_x8=False),
    ffrs.RSi16md(256, ecc_len=128, simd_x16=False, simd_x8=False),
    # ffrs.RSi16md(1024, ecc_len=128, simd_x16=False, simd_x8=False),
    # ffrs.RSi16md(4096, ecc_len=512, simd_x16=False, simd_x8=False),

    ffrs.RSi16md(4, ecc_len=2, simd_x16=False, simd_x8=False, simd_x4=False),
    ffrs.RSi16md(16, ecc_len=4, simd_x16=False, simd_x8=False, simd_x4=False),
    ffrs.RSi16md(128, ecc_len=2, simd_x16=False, simd_x8=False, simd_x4=False),
    ffrs.RSi16md(256, ecc_len=128, simd_x16=False, simd_x8=False, simd_x4=False),
    # ffrs.RSi16md(1024, ecc_len=128, simd_x16=False, simd_x8=False, simd_x4=False),
    # ffrs.RSi16md(4096, ecc_len=512, simd_x16=False, simd_x8=False, simd_x4=False),
])
class TestRSPower2(BaseTestRS):
    pass


@pytest.mark.parametrize("rs", [
    ffrs.RSi16md(16 * 13, ecc_len=8),
    ffrs.RSi16md(256 * 3, ecc_len=8),
])
class TestRSMult(BaseTestRS):
    # Limited support for multiples of powers of 2
    # TODO: partial intt or switch to forney algo
    test_repair = pytest.mark.skip(BaseTestRS.test_repair)
