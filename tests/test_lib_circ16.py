
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


class BaseTestCIRC:
    @staticmethod
    def _rsi_ecc_size(rs):
        return rs.outer_message_len * rs.inner_ecc_size * rs.outer_interleave

    @staticmethod
    def _rso_ecc_size(rs):
        return rs.inner_message_len * rs.outer_ecc_size * rs.outer_interleave

    def test_properties(self, rs):
        assert rs.ecc_len == rs.inner_block_len * rs.outer_block_len * rs.outer_interleave - rs.message_len
        assert rs.block_len == rs.message_len + rs.ecc_len

    def test_circ_encode(self, rs):
        buf = randbytes(rs.message_size)
        res = rs.encode(buf)
        assert len(res) == rs.ecc_size

        res_i = rs.rsi.encode(buf)
        assert len(res_i) == self._rsi_ecc_size(rs)
        assert res[:len(res_i)] == res_i

        res_o = rs.rso.encode(buf)
        assert len(res_o) == self._rso_ecc_size(rs)
        assert res[len(res_i):][:len(res_o)] == res_o

        res_io = rs.rsi.encode(res_o)
        assert len(res_io) == rs.outer_ecc_len * rs.inner_ecc_size * rs.outer_interleave

        # res_o does not encode the value 0x10000
        # limit to small inputs to reduce chance of failure
        if rs.inner_message_len < 64:
            assert res[-len(res_io):] == res_io

    def test_circ_repair(self, rs):
        buf = randbytes(rs.message_size)
        ecc = rs.encode(buf)
        assert len(ecc) == rs.ecc_size

        msg_orig = bytearray(buf)
        ecc_orig = bytearray(ecc)

        # Corrupt message
        for i in range((rs.outer_ecc_len - 1) * rs.inner_message_size):
            buf[i] = 0

        # Corrupt rsi corresponding to message
        for i in range(rs.inner_ecc_size * (rs.outer_ecc_len - 1)):
            # FIXME: repair fails if set `ecc[i] = 0`
            ecc[i] = 1

        # Corrupt rso
        for i in range(self._rsi_ecc_size(rs), self._rsi_ecc_size(rs) + rs.inner_ecc_size):
            ecc[i] = 0

        # Corrupt rsio corresponding to rso
        rsio_size = self._rsi_ecc_size(rs) + self._rso_ecc_size(rs)
        for i in range(rsio_size, rsio_size + rs.inner_ecc_size):
            # FIXME: repair fails if set `ecc[i] = 0`
            ecc[i] = 1

        rs.repair(buf, ecc)

        assert buf == msg_orig
        assert ecc == ecc_orig

    def test_repair_fallback(self, rs):
        buf = randbytes(rs.message_size)
        ecc = rs.encode(buf)
        assert len(ecc) == rs.ecc_size

        msg_orig = bytearray(buf)
        ecc_orig = bytearray(ecc)

        # Corrupt all rsi ecc
        for i in range(rs.inner_ecc_size * rs.outer_message_len):
            ecc[i] = 0

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


@pytest.mark.parametrize("rs", [
    ffrs.CIRC(4, 2, 4, 2),
    ffrs.CIRC(8, 2, 8, 2),
    # ffrs.CIRC(8, 2, 65536, 2),  # TODO: validate RSi16md with 64k block
    ffrs.CIRC(16, 4, 16, 4),
    ffrs.CIRC(32, 16, 4, 2),
    ffrs.CIRC(128, 8, 128, 128 // 16),
    ffrs.CIRC(512, 4, 512, 512 // 16),
    ffrs.CIRC(32768, 4, 16, 4),
    ffrs.CIRC(16, 4, 32768, 4),
    # ffrs.CIRC(1024, 4, 1024, 1024 // 16),
    # ffrs.CIRC(4096, 4, 4096, 4096 // 16),

    ffrs.CIRC(4, 2, 4, 2, simd_x16=False),
    ffrs.CIRC(8, 2, 8, 2, simd_x16=False),
    ffrs.CIRC(16, 4, 16, 4, simd_x16=False),
    ffrs.CIRC(32, 16, 4, 2, simd_x16=False),
    ffrs.CIRC(128, 8, 128, 128 // 16, simd_x16=False),
    ffrs.CIRC(512, 4, 512, 512 // 16, simd_x16=False),
    # ffrs.CIRC(1024, 4, 1024, 1024 // 16, simd_x16=False),
    # ffrs.CIRC(4096, 4, 4096, 4096 // 16, simd_x16=False),

    ffrs.CIRC(4, 2, 4, 2, simd_x16=False, simd_x8=False),
    ffrs.CIRC(8, 2, 8, 2, simd_x16=False, simd_x8=False),
    ffrs.CIRC(16, 4, 16, 4, simd_x16=False, simd_x8=False),
    ffrs.CIRC(32, 16, 4, 2, simd_x16=False, simd_x8=False),
    ffrs.CIRC(128, 8, 128, 128 // 16, simd_x16=False, simd_x8=False),
    ffrs.CIRC(512, 4, 512, 512 // 16, simd_x16=False, simd_x8=False),
    # ffrs.CIRC(1024, 4, 1024, 1024 // 16, simd_x16=False, simd_x8=False),
    # ffrs.CIRC(4096, 4, 4096, 4096 // 16, simd_x16=False, simd_x8=False),

    ffrs.CIRC(4, 2, 4, 2, simd_x16=False, simd_x8=False, simd_x4=False),
    ffrs.CIRC(8, 2, 8, 2, simd_x16=False, simd_x8=False, simd_x4=False),
    ffrs.CIRC(16, 4, 16, 4, simd_x16=False, simd_x8=False, simd_x4=False),
    ffrs.CIRC(32, 16, 4, 2, simd_x16=False, simd_x8=False, simd_x4=False),
    ffrs.CIRC(128, 8, 128, 128 // 16, simd_x16=False, simd_x8=False, simd_x4=False),
    ffrs.CIRC(512, 4, 512, 512 // 16, simd_x16=False, simd_x8=False, simd_x4=False),
    # ffrs.CIRC(1024, 4, 1024, 1024 // 16, simd_x16=False, simd_x8=False, simd_x4=False),
    # ffrs.CIRC(4096, 4, 4096, 4096 // 16, simd_x16=False, simd_x8=False, simd_x4=False),

    ffrs.CIRC(128, 2, 32, 2, 32),
    ffrs.CIRC(512, 4, 128, 8, 8),
    ffrs.CIRC(128, 8, 128, 128 // 16, 3),
    ffrs.CIRC(128, 8, 128, 128 // 16, 4),
    ffrs.CIRC(128, 8, 128, 128 // 16, 4, simd_x16=False),
    ffrs.CIRC(128, 8, 128, 128 // 16, 4, simd_x16=False, simd_x8=False),
    ffrs.CIRC(128, 8, 128, 128 // 16, 4, simd_x16=False, simd_x8=False, simd_x4=False),
])
class TestCircPower2(BaseTestCIRC):
    pass


@pytest.mark.parametrize("rs", [
    ffrs.CIRC(4*3, 2, 8*3, 2),
    ffrs.CIRC(100, 2, 100, 2),
    ffrs.CIRC(1026, 2, 128*3, 128),
    ffrs.CIRC(100, 2, 100, 2, 3),
    ffrs.CIRC(1026, 2, 128*3, 128, 7),
])
class TestCircMult(BaseTestCIRC):
    # Limited support for multiples of powers of 2
    # TODO: partial intt or switch to forney algo
    test_circ_repair = pytest.mark.skip(BaseTestCIRC.test_circ_repair)
    test_repair_fallback = pytest.mark.skip(BaseTestCIRC.test_repair_fallback)
    test_repair_no_errors = pytest.mark.skip(BaseTestCIRC.test_repair_no_errors)
