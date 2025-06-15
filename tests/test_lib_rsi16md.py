
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
        vec_size = 4
        buf = randbytes(rs.message_len * vec_size)

        res = rs.encode(buf)

        rs1 = ffrs.RSi16(rs.block_size, ecc_len=rs.ecc_len)
        for i in range(vec_size):
            ref_block = buf[i * rs.message_len:(i + 1) * rs.message_len] + bytearray(rs.ecc_len)
            ref = rs1.encode(ref_block)
            assert res[i * rs.block_size:(i + 1) * rs.block_size] == ref

    def test_encode_blocks_single(self, rs):
        vec_size = 4
        buf = randbytes(rs.message_len * vec_size)

        buf_enc = rs.encode(buf)
        buf_enc_blk = rs.encode_blocks(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.block_size * vec_size
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize('count', [3, 4, 5, 6, 7, 8, 9, 12, 16, 100])
    def test_encode_blocks_multiple(self, rs, count):
        vec_size = 4
        buf = randbytes(rs.message_len * vec_size * count)

        buf_enc = [rs.encode(buf[i * rs.message_len * vec_size:(i + 1) * rs.message_len * vec_size])
                        for i in range(count)]
        buf_enc = b"".join(buf_enc)
        buf_enc_blk = rs.encode_blocks(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.block_size * vec_size * count
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
        vec_size = 4
        buf = randbytes(rs.message_len * vec_size)

        res = rs.encode(buf)

        rs1 = ffrs.RSi16(rs.block_size, ecc_len=rs.ecc_len)
        for i in range(vec_size):
            ref_block = buf[i * rs.message_len:(i + 1) * rs.message_len] + bytearray(rs.ecc_len)
            ref = rs1.encode(ref_block)
            assert res[i * rs.ecc_len:(i + 1) * rs.ecc_len] == ref[-rs.ecc_len:]

    def test_encode_blocks_empty(self, rs):
        assert bytearray() == rs.encode_blocks(bytearray())

    def test_encode_blocks_single(self, rs):
        vec_size = 4
        buf = randbytes(rs.message_len * vec_size)

        buf_enc = rs.encode(buf)
        buf_enc_blk = rs.encode_blocks(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_len * vec_size
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize('count', [3, 4, 5, 6, 7, 8, 9, 12, 16, 100])
    def test_encode_blocks_multiple(self, rs, count):
        vec_size = 4
        buf = randbytes(rs.message_len * vec_size * count)

        buf_enc = [rs.encode(buf[i * rs.message_len * vec_size:(i + 1) * rs.message_len * vec_size])
                        for i in range(count)]
        buf_enc = b"".join(buf_enc)
        buf_enc_blk = rs.encode_blocks(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_len * vec_size * count
        assert buf_enc == buf_enc_blk

    # def test_encode_blocks_remainder(self, rs):
    #     TODO
