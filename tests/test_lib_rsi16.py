
#  test_lib_rsi16.py
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

import math
import pytest
import random

import ffrs
import ffrs.reference
import ffrs.reference.ntt
from ffrs.reference.util import to_int_list, to_bytearray, rbo_sorted

# random.seed(42)

gfref = ffrs.reference.GF(65537, 1, 3)
GF = ffrs.GFi16(65537, 3)


@pytest.mark.parametrize('rs', [
    ffrs.RSi16(8, ecc_len=2*2),
    ffrs.RSi16(32, ecc_len=4*2),
    ffrs.RSi16(256, ecc_len=8*2),
    ffrs.RSi16(256, ecc_len=9*2),
    ffrs.RSi16(256, ecc_len=10*2),
    ffrs.RSi16(512, ecc_len=4*2),
])
class TestRS:
    def test_encode(self, rs):
        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        buf = [random.randrange(0, 2**16) for _ in range(size_u16 - ecc_u16)] + [0] * ecc_u16

        w = gfref(rs.roots_of_unity[round(math.log2(size_u16))])

        ref = to_int_list(ffrs.reference.ntt.ntt(gfref, gfref(w), rbo_sorted(buf)))
        res_i = to_int_list(ffrs.reference.ntt.intt(gfref, gfref(w), ref))

        res = to_bytearray(buf, 2)
        rs.encode(res)
        res = to_int_list(res, 2)

        res2 = rs.gf.ntt16(to_bytearray(buf, 2), w)
        res2 = to_int_list(res2, 2)

        assert rbo_sorted(res_i) == buf
        assert ref[:ecc_u16] == res[-ecc_u16:]

    def test_encode_decode(self, rs):
        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        orig = [random.randrange(0, 2**16) for _ in range(size_u16 - ecc_u16)] + [0] * ecc_u16

        msg_enc = to_bytearray(orig, 2)
        rs.encode(msg_enc)

        msg_enc_err = bytearray(msg_enc)
        assert msg_enc == msg_enc_err
        msg_enc_err = to_int_list(msg_enc_err, 2)

        err_vec = [0] * size_u16

        error_positions = dict()
        while len(error_positions) < ecc_u16//2:
            # Do not corrupt codeword for now
            i = random.randrange(size_u16 - ecc_u16)
            if i in error_positions:
                continue
            error_positions[i] = random.randrange(2**16)
            msg_enc_err[i] = rs.gf.add(msg_enc_err[i], error_positions[i])
            err_vec[i] = rs.gf.add(err_vec[i], error_positions[i])

        assert msg_enc_err != msg_enc

        msg_enc_err_dec = to_bytearray(msg_enc_err, 2)
        rs.decode(msg_enc_err_dec)

        assert msg_enc == msg_enc_err_dec

    def test_encode_blocks(self, rs):
        blocks = 3

        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        msg_u16 = rs.message_len // 2
        buf = [random.randrange(1, 2**16) for _ in range(msg_u16 * blocks)]

        w = gfref(rs.roots_of_unity[round(math.log2(size_u16))])

        buf_ntt = [to_int_list(ffrs.reference.ntt.ntt(gfref, gfref(w), rbo_sorted(buf[msg_u16 * blk:msg_u16 * (blk + 1)] + [0] * ecc_u16))) for blk in range(blocks)]
        res_i = [to_int_list(ffrs.reference.ntt.intt(gfref, gfref(w), rbo_sorted(buf_ntt[blk]))) for blk in range(blocks)]

        res = rs.encode_blocks(to_bytearray(buf, 2))
        res = to_int_list(res, 2)

        for i in range(blocks):
            start = (i + 1) * size_u16 - ecc_u16
            assert buf_ntt[i][:ecc_u16] == res[start:start + ecc_u16]

