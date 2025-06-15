
#  test_lib_rsi16md.py
#
#  Copyright 2025 Gabriel Machado
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
import ffrs.reference
import ffrs.reference.ntt
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted, randbytes

# random.seed(42)

ref_gf = ffrs.reference.GF(65537, 1, 3)

ref_ntt_ct = lambda w, buf: to_int_list(ffrs.reference.ntt.ntt(ref_gf, ref_gf(w), rbo_sorted(buf)))
ref_intt_ct = lambda w, buf: rbo_sorted(to_int_list(ffrs.reference.ntt.intt(ref_gf, ref_gf(w), buf)))
ref_ntt_gs = lambda w, buf: rbo_sorted(to_int_list(ffrs.reference.ntt.ntt(ref_gf, ref_gf(w), buf)))
ref_intt_gs = lambda w, buf: to_int_list(ffrs.reference.ntt.intt(ref_gf, ref_gf(w), rbo_sorted(buf)))

ref_ntt, ref_intt = ref_ntt_ct, ref_intt_ct
# ref_ntt, ref_intt = ref_ntt_gs, ref_intt_gs


@pytest.mark.parametrize('rs', [
    ffrs.RSi16md(4*2, ecc_len=2*2, inline=True),
    ffrs.RSi16md(16*2, ecc_len=4*2, inline=True),
    ffrs.RSi16md(128*2, ecc_len=2*2, inline=True),
    ffrs.RSi16md(256*2, ecc_len=128*2, inline=True),
    ffrs.RSi16md(1024, ecc_len=128, inline=True),
])
class TestRSInline:
    def test_encode(self, rs):
        buf = randbytes(rs.message_len)

        res = rs.encode(buf)

        # Compare inline/external
        rs2 = ffrs.RSi16md(rs.block_size, rs.message_len, inline=False)
        res_ecc = rs2.encode(buf)
        assert res[-rs.ecc_len:] == res_ecc

        w = ref_gf(rs.root)
        buf_ntt = ref_ntt(w, to_int_list(buf + bytearray(rs.ecc_len), 2))

        ecc_mix = rs._mix_ecc(w, buf_ntt[:rs.ecc_len//2])

        assert to_bytearray(ecc_mix, 2) == res[-rs.ecc_len:]
        assert res[:rs.message_len] == buf

    def test_encode_decode(self, rs):
        size_u16 = rs.block_size // 2
        ecc_u16 = rs.ecc_len // 2
        orig = [random.randrange(0, 2**16) for _ in range(size_u16 - ecc_u16)]

        msg_enc = to_int_list(rs.encode(to_bytearray(orig, 2)), 2)

        msg_enc_err = list(msg_enc)

        assert rs.find_errors(to_bytearray(msg_enc, 2)) == {}

        error_positions = dict()

        # Always add one error to the codeword
        if len(error_positions) < max(ecc_u16//2, 1):
            i = random.randrange(size_u16 - ecc_u16, size_u16)
            error_positions[i] = random.randrange(2**16)
            msg_enc_err[i] = rs.gf.add(msg_enc_err[i], error_positions[i])

        while len(error_positions) < ecc_u16//2:
            i = random.randrange(size_u16)
            if i in error_positions:
                continue
            error_positions[i] = random.randrange(2**16)
            msg_enc_err[i] = rs.gf.add(msg_enc_err[i], error_positions[i])

        w = ref_gf(rs.root)
        print()
        print("errs: ", error_positions)

        assert msg_enc_err != msg_enc

        msg_enc_err_dec = rs.decode(to_bytearray(msg_enc_err, 2))

        assert msg_enc[:-ecc_u16] == to_int_list(msg_enc_err_dec, 2)
        assert error_positions == rs.find_errors(to_bytearray(msg_enc_err, 2))

    def test_encode_blocks_single(self, rs):
        buf = randbytes(rs.message_len)

        buf_enc = rs.encode(buf)
        buf_enc_blk = rs.encode_blocks(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.block_size
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize('count', [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 100])
    def test_encode_blocks_multiple(self, rs, count):
        buf = randbytes(rs.message_len * count)

        buf_enc = [rs.encode(buf[i * rs.message_len:(i + 1) * rs.message_len])
                        for i in range(count)]
        buf_enc = b"".join(buf_enc)
        buf_enc_blk = rs.encode_blocks(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.block_size * count
        assert buf_enc == buf_enc_blk

    # def test_encode_blocks_remainder(self, rs):
    #     TODO


@pytest.mark.parametrize('rs', [
    ffrs.RSi16md(4*2, ecc_len=2*2),
    ffrs.RSi16md(16*2, ecc_len=4*2),
    ffrs.RSi16md(128*2, ecc_len=2*2),
    ffrs.RSi16md(256*2, ecc_len=128*2),
    ffrs.RSi16md(1024, ecc_len=128),
])
class TestRSExternal:
    def test_encode(self, rs):
        buf = randbytes(rs.message_len)

        res = rs.encode(buf)

        w = ref_gf(rs.root)
        buf_ntt = ref_ntt(w, to_int_list(buf + bytearray(rs.ecc_len), 2))

        ecc_mix = rs._mix_ecc(w, buf_ntt[:rs.ecc_len//2])

        assert to_bytearray(ecc_mix, 2) == res

    def test_encode_blocks_empty(self, rs):
        assert bytearray() == rs.encode_blocks(bytearray())

    def test_encode_blocks_single(self, rs):
        buf = randbytes(rs.message_len)

        buf_enc = rs.encode(buf)
        buf_enc_blk = rs.encode_blocks(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_len
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize('count', [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 100])
    def test_encode_blocks_multiple(self, rs, count):
        buf = randbytes(rs.message_len * count)

        buf_enc = [rs.encode(buf[i * rs.message_len:(i + 1) * rs.message_len])
                        for i in range(count)]
        buf_enc = b"".join(buf_enc)
        buf_enc_blk = rs.encode_blocks(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_len * count
        assert buf_enc == buf_enc_blk

    # def test_encode_blocks_remainder(self, rs):
    #     TODO
