
#  test_lib_circ16.py
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

import ffrs

from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted, randbytes


@pytest.mark.parametrize("rs", [
    ffrs.CIRC(4, 2, 4, 2),
    ffrs.CIRC(32, 2, 32, 4),
    ffrs.CIRC(128, 4, 128, 128 // 16),
    ffrs.CIRC(512, 4, 512, 512 // 16),
    ffrs.CIRC(1024, 4, 1024, 1024 // 16),
    # ffrs.CIRC(4096, 4, 4096, 4096 // 16),
])
class TestCIRC:
    def test_circ_encode(self, rs):
        buf = randbytes(rs.message_size)
        res = rs.encode(buf)
        assert len(res) == rs.ecc_size
        assert len(res) // 2 == rs.rsi.block_len * rs.rso.block_len - rs.message_len

        res_i = rs.rsi.encode_blocks(buf)
        assert len(res_i) == rs.rso.message_len * rs.rsi.ecc_size
        assert res[:len(res_i)] == res_i

        res_o = rs.rso.encode_chunk(buf)
        assert len(res_o) == rs.rsi.message_len * rs.rso.ecc_size
        assert res[len(res_i):][:len(res_o)] == res_o

    def test_circ_repair(self, rs):
        buf = randbytes(rs.message_size)
        ecc = rs.encode(buf)
        assert len(ecc) == rs.ecc_size

        msg_orig = bytearray(buf)
        ecc_orig = bytearray(ecc)

        # TODO: repair errors in RSO ECC
        grace = 2

        for i in range((rs.rso.ecc_len - grace) * rs.rsi.message_size):
            buf[i] = 0

        rs.repair(buf, ecc)

        assert buf == msg_orig
        assert ecc == ecc_orig

    def test_repair_no_errors(self, rs):
        buf = bytearray(randbytes(rs.message_size))
        ecc = rs.encode(buf)

        msg_orig = bytearray(buf)
        ecc_orig = bytearray(ecc)

        rs.repair(buf, ecc)

        assert buf == msg_orig
        assert ecc == ecc_orig
