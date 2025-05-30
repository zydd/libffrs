
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


from common import to_int_list, to_bytearray

random.seed(42)

GF65537 = ffrs.reference.GF(65537, 1, 3)
GF = ffrs.GFi32(65537, 3)


@pytest.mark.parametrize('rs', [
    ffrs.RSi16(16, ecc_len=4),
    ffrs.RSi16(256, ecc_len=200),
    ffrs.RSi16(512, ecc_len=20),
])
class TestRS:
    def test_encode(self, rs):
        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        buf = [random.randrange(0, 2**16) for _ in range(size_u16 - ecc_u16)] + [0] * ecc_u16

        w = GF65537(rs.roots_of_unity[round(math.log2(size_u16))])

        ref = to_int_list(ffrs.reference.ntt.ntt(GF65537, GF65537(w), buf))
        res_i = to_int_list(ffrs.reference.ntt.intt(GF65537, GF65537(w), ref))

        res = to_bytearray(buf, 2)
        rs.encode(res)
        res = to_int_list(res, 2)

        res2 = rs.gf.ntt16(to_bytearray(buf, 2), w)
        res2 = to_int_list(res2, 2)

        assert res_i == buf
        assert ref[:ecc_u16] == res[-ecc_u16:]

    def test_encode_blocks(self, rs):
        blocks = 3

        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        msg_u16 = rs.message_len // 2
        buf = [random.randrange(1, 2**16) for _ in range(msg_u16 * blocks)]

        w = GF65537(rs.roots_of_unity[round(math.log2(size_u16))])

        buf_ntt = [to_int_list(ffrs.reference.ntt.ntt(GF65537, GF65537(w), buf[msg_u16 * blk:msg_u16 * (blk + 1)] + [0] * ecc_u16)) for blk in range(blocks)]
        res_i = [to_int_list(ffrs.reference.ntt.intt(GF65537, GF65537(w), buf_ntt[blk])) for blk in range(blocks)]
        res_i = [item for sublist in res_i for item in sublist[:-ecc_u16]]

        res = rs.encode_blocks(to_bytearray(buf, 2))
        res = to_int_list(res, 2)

        assert res_i == buf
        for i in range(blocks):
            start = (i + 1) * size_u16 - ecc_u16
            assert buf_ntt[i][:ecc_u16] == res[start:start + ecc_u16]
