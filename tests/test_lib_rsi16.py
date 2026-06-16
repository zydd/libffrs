#  test_lib_rsi16.py
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


def add_errors(msg, ecc, count):
    error_positions = {}

    # # Force 1 error in ecc
    # i = random.randrange(len(ecc))
    # error_positions[len(msg) + i] = random.randrange(2**16)
    # ecc[i] ^= error_positions[len(msg) + i]

    while len(error_positions) < count:
        i = random.randrange((len(msg) + len(ecc)))
        if i in error_positions:
            continue
        error_positions[i] = random.randrange(1, 256)
        if i < len(msg):
            msg[i] = msg[i] ^ error_positions[i]
        else:
            ecc[i - len(msg)] = ecc[i - len(msg)] ^ error_positions[i]

    return msg, ecc, error_positions


class BaseTestRS:
    def test_encode(self, rs: ffrs.RSi16):
        buf = randbytes(rs.message_size)

        res = rs.encode(buf)

        w = ref_gf(rs.root)
        buf_ntt = ref_ntt(w, to_int_list(buf + bytearray(rs.ntt_size - len(buf))))

        ecc_mix = ref_rs.mix_ecc(rs, ref_gf, buf_ntt[: rs.ecc_len])

        assert to_bytearray(ecc_mix) == res

    @pytest.mark.parametrize("grace", [0, 1, 2, 3, 4])
    def test_repair_unknown_locations(self, rs: ffrs.RSi16, grace):
        msg_orig = [random.randrange(0, 2**16) for _ in range(rs.message_len)]
        ecc_orig = to_int_list(rs.encode(to_bytearray(msg_orig)))
        msg_err = list(msg_orig)
        ecc_err = list(ecc_orig)

        msg_err, ecc_err, _error_positions = add_errors(msg_err, ecc_err, max(rs.ecc_len // 2 - grace, 1))

        assert msg_err != msg_orig or ecc_err != ecc_orig

        msg_err = to_bytearray(msg_err)
        ecc_err = to_bytearray(ecc_err)
        rs.repair(msg_err, ecc_err)
        msg_err = to_int_list(msg_err)
        ecc_err = to_int_list(ecc_err)

        assert msg_err + ecc_err == msg_orig + ecc_orig

    @pytest.mark.parametrize("grace", [0, 1, 2, 3, 4])
    def test_repair(self, rs: ffrs.RSi16, grace):
        msg_orig = [random.randrange(0, 2**16) for _ in range(rs.message_len)]
        ecc_orig = to_int_list(rs.encode(to_bytearray(msg_orig)))
        msg_err = list(msg_orig)
        ecc_err = list(ecc_orig)

        msg_err, ecc_err, error_positions = add_errors(msg_err, ecc_err, max(rs.ecc_len - grace, 1))

        assert msg_err + ecc_err != msg_orig + ecc_orig

        msg_err = to_bytearray(msg_err)
        ecc_err = to_bytearray(ecc_err)
        rs.repair(msg_err, ecc_err, list(error_positions.keys()))
        msg_err = to_int_list(msg_err)
        ecc_err = to_int_list(ecc_err)

        assert msg_err + ecc_err == msg_orig + ecc_orig

    def test_repair_no_errors(self, rs: ffrs.RSi16):
        orig = [random.randrange(0, 2**16) for _ in range(rs.block_len - rs.ecc_len)]
        msg = to_bytearray(orig)
        ecc = rs.encode(msg)
        msg_err = bytearray(msg)
        ecc_err = bytearray(ecc)

        rs.repair(msg_err, ecc_err, [])

        assert msg_err == msg
        assert ecc_err == ecc

    def test_repair_unknown_no_errors(self, rs: ffrs.RSi16):
        orig = [random.randrange(0, 2**16) for _ in range(rs.block_len - rs.ecc_len)]
        msg = to_bytearray(orig)
        ecc = rs.encode(msg)
        msg_err = bytearray(msg)
        ecc_err = bytearray(ecc)

        rs.repair(msg_err, ecc_err)

        assert msg_err == msg
        assert ecc_err == ecc

    def test_encode_blocks_empty(self, rs: ffrs.RSi16):
        assert bytearray() == rs.encode(bytearray())

    def test_encode_blocks_single(self, rs: ffrs.RSi16):
        buf = randbytes(rs.message_size)

        buf_enc = rs.encode(buf)
        buf_enc_blk = rs.encode(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_size
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize("count", [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 100])
    def test_encode_blocks_multiple(self, rs: ffrs.RSi16, count):
        buf = randbytes(rs.message_size * count)

        buf_enc = [rs.encode(buf[i * rs.message_size : (i + 1) * rs.message_size]) for i in range(count)]
        buf_enc = b"".join(buf_enc)
        buf_enc_blk = rs.encode(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_size * count
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize("extra", [2, 4, 6, 8, 10, 12, 14, 128, 1500])
    @pytest.mark.skip("not supported for now")
    def test_encode_blocks_remainder(self, rs: ffrs.RSi16, extra):
        count = 16 + extra // rs.message_size
        extra = (extra % rs.message_size) or 2
        buf = randbytes(rs.message_size * count + extra)

        buf_enc = [rs.encode(buf[i * rs.message_size : (i + 1) * rs.message_size]) for i in range(count)]

        buf_enc += [rs.encode(buf[-extra:] + bytearray(rs.message_size - extra))]
        buf_enc = b"".join(buf_enc)
        buf_enc_blk = rs.encode(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_size * (count + 1)
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize("interleave", list(range(2, 15)) + [32, 48, 100, 256])
    def test_encode_interleaved(self, rs: ffrs.RSi16, interleave):
        assert rs.interleave == 1
        rsi = ffrs.RSi16(
            rs.block_len,
            rs.ecc_len,
            interleave=interleave,
            simd_x4=rs.simd_x4,
            simd_x8=rs.simd_x8,
            simd_x16=rs.simd_x16,
        )

        data = list(range(interleave * rs.message_len))
        buf = to_bytearray(data)

        interleaved_chunk = rsi.encode(buf)

        for i in range(interleave):
            buf_enc = rs.encode(to_bytearray(data[i::interleave]))

            # Sequential ecc
            # assert interleaved_chunk[i * rs.ecc_size:(i + 1) * rs.ecc_size] == buf_enc

            # Interleaved ecc
            assert to_bytearray(to_int_list(interleaved_chunk)[i::interleave]) == buf_enc

    @pytest.mark.parametrize("interleave", list(range(1, 16)) + [32, 48, 50])
    @pytest.mark.parametrize("grace", range(4))
    def test_repair_interleaved_unknown_unaligned(self, rs: ffrs.RSi16, interleave, grace):
        assert rs.interleave == 1
        rsi = ffrs.RSi16(
            rs.block_len,
            rs.ecc_len,
            interleave=interleave,
            simd_x4=rs.simd_x4,
            simd_x8=rs.simd_x8,
            simd_x16=rs.simd_x16,
        )

        msg_orig = list(range(rsi.message_len))
        msg_buf_orig = to_bytearray(msg_orig)

        ecc_buf_orig = rsi.encode(msg_buf_orig)
        ecc_orig = to_int_list(ecc_buf_orig)

        msg_err = list(msg_orig)
        ecc_err = list(ecc_orig)

        errors = []
        for i in range(interleave):
            msg_err[i::interleave], ecc_err[i::interleave], col_errs = add_errors(
                msg_err[i::interleave], ecc_err[i::interleave], max(rs.ecc_len // 2 - grace, 1)
            )
            errors.append(col_errs)

        msg_err0 = msg_err[::interleave]
        ecc_err0 = ecc_err[::interleave]
        rs.repair(to_bytearray(msg_err0), to_bytearray(ecc_err0))
        msg_buf_err = to_bytearray(msg_err)
        ecc_buf_err = to_bytearray(ecc_err)

        assert msg_buf_err + ecc_buf_err != msg_buf_orig + ecc_buf_orig
        rsi.repair(msg_buf_err, ecc_buf_err)
        assert msg_buf_err + ecc_buf_err == msg_buf_orig + ecc_buf_orig

    @pytest.mark.parametrize("interleave", list(range(1, 16)) + [32, 48, 50])
    @pytest.mark.parametrize("grace", range(4))
    def test_repair_interleaved_unknown_aligned(self, rs: ffrs.RSi16, interleave, grace):
        assert rs.interleave == 1
        rsi = ffrs.RSi16(
            rs.block_len,
            rs.ecc_len,
            interleave=interleave,
            simd_x4=rs.simd_x4,
            simd_x8=rs.simd_x8,
            simd_x16=rs.simd_x16,
        )

        msg_orig = list(range(rsi.message_len))
        msg_buf_orig = to_bytearray(msg_orig)

        ecc_buf_orig = rsi.encode(msg_buf_orig)
        ecc_orig = to_int_list(ecc_buf_orig)

        msg_err = list(msg_orig)
        ecc_err = list(ecc_orig)

        # Aligned errors
        msg_err[::interleave], ecc_err[::interleave], errors = add_errors(
            msg_err[::interleave], ecc_err[::interleave], max(rs.ecc_len // 2 - grace, 1)
        )
        for i in range(1, interleave):
            for pos, _err in errors.items():
                if pos < rs.message_len:
                    msg_err[i + pos * interleave] ^= random.randint(1, 65536)
                else:
                    ecc_err[i + (pos - rs.message_len) * interleave] ^= random.randint(1, 65536)

        msg_err0 = msg_err[::interleave]
        ecc_err0 = ecc_err[::interleave]
        rs.repair(to_bytearray(msg_err0), to_bytearray(ecc_err0))

        msg_buf_err = to_bytearray(msg_err)
        ecc_buf_err = to_bytearray(ecc_err)

        assert msg_buf_err + ecc_buf_err != msg_buf_orig + ecc_buf_orig
        rsi.repair(msg_buf_err, ecc_buf_err)
        assert msg_buf_err + ecc_buf_err == msg_buf_orig + ecc_buf_orig

    @pytest.mark.parametrize("interleave", list(range(1, 16)) + [32, 33, 48, 50])
    @pytest.mark.parametrize("grace", range(4))
    def test_repair_interleaved_known_aligned(self, rs: ffrs.RSi16, interleave, grace):
        assert rs.interleave == 1
        rsi = ffrs.RSi16(
            rs.block_len,
            rs.ecc_len,
            interleave=interleave,
            simd_x4=rs.simd_x4,
            simd_x8=rs.simd_x8,
            simd_x16=rs.simd_x16,
        )

        msg_orig = list(range(rsi.message_len))
        msg_buf_orig = to_bytearray(msg_orig)

        ecc_buf_orig = rsi.encode(msg_buf_orig)
        ecc_orig = to_int_list(ecc_buf_orig)

        msg_err = list(msg_orig)
        ecc_err = list(ecc_orig)

        # Aligned errors
        msg_err[::interleave], ecc_err[::interleave], errors = add_errors(
            msg_err[::interleave], ecc_err[::interleave], max(rs.ecc_len - grace, 1)
        )
        for i in range(1, interleave):
            for pos, _err in errors.items():
                if pos < rs.message_len:
                    msg_err[i + pos * interleave] ^= random.randint(1, 65536)
                else:
                    ecc_err[i + (pos - rs.message_len) * interleave] ^= random.randint(1, 65536)

        msg_err0 = to_bytearray(msg_err[::interleave])
        ecc_err0 = to_bytearray(ecc_err[::interleave])
        rs.repair(msg_err0, ecc_err0, errors.keys())
        assert to_int_list(msg_err0) == msg_orig[::interleave]
        assert to_int_list(ecc_err0) == ecc_orig[::interleave]

        msg_buf_err = to_bytearray(msg_err)
        ecc_buf_err = to_bytearray(ecc_err)

        assert msg_buf_err + ecc_buf_err != msg_buf_orig + ecc_buf_orig
        rsi.repair(msg_buf_err, ecc_buf_err, errors.keys())
        assert msg_buf_err + ecc_buf_err == msg_buf_orig + ecc_buf_orig


@pytest.mark.parametrize(
    "rs",
    [
        ffrs.RSi16(4, ecc_len=2),
        ffrs.RSi16(8, ecc_len=2),
        ffrs.RSi16(16, ecc_len=4),
        ffrs.RSi16(16, ecc_len=8),
        ffrs.RSi16(128, ecc_len=2),
        ffrs.RSi16(256, ecc_len=32),
        # ffrs.RSi16(1024, ecc_len=128),
        # ffrs.RSi16(4096, ecc_len=512),
        ffrs.RSi16(4, ecc_len=2, simd_x16=False),
        ffrs.RSi16(16, ecc_len=4, simd_x16=False),
        ffrs.RSi16(128, ecc_len=2, simd_x16=False),
        ffrs.RSi16(256, ecc_len=32, simd_x16=False),
        # ffrs.RSi16(1024, ecc_len=128, simd_x16=False),
        # ffrs.RSi16(4096, ecc_len=512, simd_x16=False),
        ffrs.RSi16(4, ecc_len=2, simd_x16=False, simd_x8=False),
        ffrs.RSi16(16, ecc_len=4, simd_x16=False, simd_x8=False),
        ffrs.RSi16(128, ecc_len=2, simd_x16=False, simd_x8=False),
        ffrs.RSi16(256, ecc_len=32, simd_x16=False, simd_x8=False),
        # ffrs.RSi16(1024, ecc_len=128, simd_x16=False, simd_x8=False),
        # ffrs.RSi16(4096, ecc_len=512, simd_x16=False, simd_x8=False),
        ffrs.RSi16(4, ecc_len=2, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.RSi16(16, ecc_len=4, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.RSi16(128, ecc_len=2, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.RSi16(256, ecc_len=32, simd_x16=False, simd_x8=False, simd_x4=False),
        # ffrs.RSi16(1024, ecc_len=128, simd_x16=False, simd_x8=False, simd_x4=False),
        # ffrs.RSi16(4096, ecc_len=512, simd_x16=False, simd_x8=False, simd_x4=False),
    ],
)
class TestRSPower2(BaseTestRS):
    pass


@pytest.mark.parametrize(
    "rs",
    [
        ffrs.RSi16(16 * 13, ecc_len=8),
        ffrs.RSi16(256 * 3, ecc_len=8),
    ],
)
class TestRSMult(BaseTestRS):
    pass


@pytest.mark.parametrize(
    "rs",
    [
        ffrs.RSi16(4096, ecc_len=8),
        ffrs.RSi16(4096, ecc_len=256),
        ffrs.RSi16(2048, ecc_len=128),
        ffrs.RSi16(4096, ecc_len=2048),
        ffrs.RSi16(32768, ecc_len=256),
        ffrs.RSi16(65536, ecc_len=256),
    ],
)
class TestRSLarge:
    def test_repair_no_errors(self, rs: ffrs.RSi16):
        msg = randbytes(rs.message_size)
        ecc = rs.encode(msg)
        msg_err = bytearray(msg)
        ecc_err = bytearray(ecc)

        rs.repair(msg_err, ecc_err, [])

        assert msg_err == msg
        assert ecc_err == ecc

    def test_repair_unknown_no_errors(self, rs: ffrs.RSi16):
        msg = randbytes(rs.message_size)
        ecc = rs.encode(msg)
        msg_err = bytearray(msg)
        ecc_err = bytearray(ecc)

        rs.repair(msg_err, ecc_err)

        assert msg_err == msg
        assert ecc_err == ecc

    @pytest.mark.parametrize("grace", [2, 3, 4, 5])
    def test_repair_unknown_locations(self, rs: ffrs.RSi16, grace):
        msg_orig = randbytes(rs.message_size)
        ecc_orig = rs.encode(msg_orig)

        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        msg_err, ecc_err, _error_positions = add_errors(msg_err, ecc_err, max(rs.ecc_len // 2 - grace, 1))

        # Corrupt first and last bytes of the block
        msg_err[-1] ^= 0xFF
        ecc_err[-1] ^= 0xFF

        assert msg_err != msg_orig or ecc_err != ecc_orig

        rs.repair(msg_err, ecc_err)

        assert msg_err == msg_orig
        assert ecc_err == ecc_orig

    @pytest.mark.parametrize("grace", range(4))
    def test_repair_known_locations(self, rs: ffrs.RSi16, grace):
        msg_orig = randbytes(rs.message_size)
        ecc_orig = rs.encode(msg_orig)

        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        msg_err, ecc_err, error_positions = add_errors(msg_err, ecc_err, max(rs.ecc_len - grace, 1))

        assert msg_err != msg_orig or ecc_err != ecc_orig

        rs.repair(msg_err, ecc_err, set(map(lambda x: x // 2, error_positions.keys())))

        assert msg_err == msg_orig
        assert ecc_err == ecc_orig
