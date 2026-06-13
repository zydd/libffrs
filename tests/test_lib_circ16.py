#
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

import logging
import random
import sys

import pytest

import ffrs
import ffrs.par

from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted, randbytes


class BaseTestCIRC:
    @staticmethod
    def _rsi_ecc_size(rs: ffrs.CIRC16):
        return rs.outer_message_len * rs.inner_ecc_size * rs.interleave

    @staticmethod
    def _rso_ecc_size(rs: ffrs.CIRC16):
        return rs.inner_message_size * rs.outer_ecc_len * rs.interleave

    @staticmethod
    def add_errors_stride(msg, ecc, stride, max=float("inf")):
        pos = 0
        count = 0
        while pos < len(msg) and count < max:
            msg[pos] ^= 0xFF
            pos += stride
            count += 1

        pos -= len(msg)
        while pos < len(ecc) and count < max:
            ecc[pos] ^= 0xFF
            pos += stride
            count += 1

        return count

    def test_properties(self, rs: ffrs.CIRC16):
        assert rs.ecc_len == rs.inner_block_len * rs.outer_block_len * rs.interleave - rs.message_len
        assert rs.block_len == rs.message_len + rs.ecc_len

    def test_circ_encode(self, rs: ffrs.CIRC16):
        buf = randbytes(rs.message_size)
        res = rs.encode(buf)
        assert len(res) == rs.ecc_size

        res_o = rs.rso.encode(buf)
        assert len(res_o) == self._rso_ecc_size(rs)
        assert res[: len(res_o)] == res_o

        res_i = rs.rsi.encode(buf)
        assert len(res_i) == self._rsi_ecc_size(rs)
        assert res[len(res_o) :][: len(res_i)] == res_i

        res_io = rs.rsi.encode(res_o)
        assert len(res_io) == rs.outer_ecc_len * rs.inner_ecc_size * rs.interleave

        # res_o does not encode the value 0x10000
        # limit to small inputs to reduce chance of failure
        if rs.inner_message_len < 64:
            assert res[-len(res_io) :] == res_io

    def add_errors_start(self, rs, buf, ecc):
        # Corrupt message
        for i in range((rs.outer_ecc_len - 1) * rs.inner_message_size):
            buf[i] ^= random.randint(0, 255)

        # Corrupt rso
        for i in range(rs.inner_ecc_size):
            ecc[i] ^= random.randint(0, 255)

        # Corrupt rsi corresponding to message
        for i in range(rs.inner_ecc_size * (rs.outer_ecc_len - 1)):
            ecc[self._rso_ecc_size(rs) + i] ^= random.randint(0, 255)

        # Corrupt rsio corresponding to rso
        rsio_size = self._rsi_ecc_size(rs) + self._rso_ecc_size(rs)
        for i in range(rsio_size, rsio_size + rs.inner_ecc_size):
            ecc[i] ^= random.randint(0, 255)

    def test_circ_repair(self, rs: ffrs.CIRC16):
        buf = randbytes(rs.message_size)
        ecc = rs.encode(buf)
        assert len(ecc) == rs.ecc_size

        msg_orig = bytearray(buf)
        ecc_orig = bytearray(ecc)

        self.add_errors_start(rs, buf, ecc)

        # buf_err = bytearray(buf)
        # ecc_err = bytearray(ecc)
        # rs2 = ffrs.CIRC16(rs.inner_block_len, rs.inner_ecc_len, rs.outer_block_len, rs.outer_ecc_len, rs.interleave, simd_x16=False, simd_x8=False, simd_x4=False)
        # rs2.repair(buf_err, ecc_err)

        rs.repair(buf, ecc)

        assert buf == msg_orig
        assert ecc == ecc_orig

    @pytest.mark.parametrize("count", range(1, 10))
    def test_circ_repair_multiple(self, rs: ffrs.CIRC16, count):
        buf = randbytes(rs.message_size)
        ecc = rs.encode(buf)
        assert len(ecc) == rs.ecc_size

        msg_orig = bytearray(buf)
        ecc_orig = bytearray(ecc)

        self.add_errors_start(rs, buf, ecc)

        rs.repair(buf, ecc)

        assert buf == msg_orig
        assert ecc == ecc_orig

    @pytest.mark.skip
    def test_circ_repair_zeroes(self, rs: ffrs.CIRC16):
        buf = randbytes(rs.message_size)
        ecc = rs.encode(buf)
        assert len(ecc) == rs.ecc_size

        msg_orig = bytearray(buf)
        ecc_orig = bytearray(ecc)

        # Corrupt message
        for i in range((rs.outer_ecc_len - 1) * rs.inner_message_size):
            buf[i] = 0

        # Corrupt rso
        for i in range(rs.inner_ecc_size):
            ecc[i] = 0

        # Corrupt rsi corresponding to message
        for i in range(rs.inner_ecc_size * (rs.outer_ecc_len - 1)):
            # FIXME: repair fails if set `ecc[i] = 0`
            # Likely because these are not marked as errors since message is also `0`
            # Then NTT repair does not work properly
            # TODO: see if this can be used in all cases:
            # https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders#Error_and_erasure_correction
            ecc[self._rso_ecc_size(rs) + i] = 0

        # Corrupt rsio corresponding to rso
        rsio_size = self._rsi_ecc_size(rs) + self._rso_ecc_size(rs)
        for i in range(rsio_size, rsio_size + rs.inner_ecc_size):
            # FIXME: repair fails if set `ecc[i] = 0`
            ecc[i] = 0

        rs.repair(buf, ecc)

        assert buf == msg_orig
        assert ecc == ecc_orig

    def test_repair_fallback(self, rs: ffrs.CIRC16):
        buf = randbytes(rs.message_size)
        ecc = rs.encode(buf)
        assert len(ecc) == rs.ecc_size

        msg_orig = bytearray(buf)
        ecc_orig = bytearray(ecc)

        # Corrupt all rsi ecc
        o = self._rso_ecc_size(rs)
        i = self._rsi_ecc_size(rs)
        ecc[o : o + i] = bytearray(self._rsi_ecc_size(rs))

        err_count = self.add_errors_stride(
            buf, [], rs.inner_message_size + 2, rs.inner_message_len * rs.outer_ecc_len // 2
        )

        rs.repair(buf, ecc)

        assert buf == msg_orig
        assert ecc == ecc_orig

    def test_repair_no_errors(self, rs: ffrs.CIRC16):
        buf = bytearray(randbytes(rs.message_size))
        ecc = rs.encode(buf)

        msg_orig = bytearray(buf)
        ecc_orig = bytearray(ecc)

        rs.repair(buf, ecc)

        assert buf == msg_orig
        assert ecc == ecc_orig


@pytest.mark.parametrize(
    "rs",
    [
        ffrs.CIRC16(4, 2, 4, 2),
        ffrs.CIRC16(8, 2, 8, 2),
        # ffrs.CIRC(8, 2, 65536, 2),  # TODO: validate RSi16 with 64k block
        ffrs.CIRC16(16, 4, 16, 4),
        ffrs.CIRC16(32, 16, 4, 2),
        ffrs.CIRC16(128, 8, 128, 128 // 16),
        ffrs.CIRC16(512, 4, 512, 512 // 16),
        ffrs.CIRC16(32768, 4, 16, 4),
        ffrs.CIRC16(16, 4, 32768, 4),
        # ffrs.CIRC(1024, 4, 1024, 1024 // 16),
        # ffrs.CIRC(4096, 4, 4096, 4096 // 16),
        ffrs.CIRC16(4, 2, 4, 2, simd_x16=False),
        ffrs.CIRC16(8, 2, 8, 2, simd_x16=False),
        ffrs.CIRC16(16, 4, 16, 4, simd_x16=False),
        ffrs.CIRC16(32, 16, 4, 2, simd_x16=False),
        ffrs.CIRC16(128, 8, 128, 128 // 16, simd_x16=False),
        ffrs.CIRC16(512, 4, 512, 512 // 16, simd_x16=False),
        # ffrs.CIRC(1024, 4, 1024, 1024 // 16, simd_x16=False),
        # ffrs.CIRC(4096, 4, 4096, 4096 // 16, simd_x16=False),
        ffrs.CIRC16(4, 2, 4, 2, simd_x16=False, simd_x8=False),
        ffrs.CIRC16(8, 2, 8, 2, simd_x16=False, simd_x8=False),
        ffrs.CIRC16(16, 4, 16, 4, simd_x16=False, simd_x8=False),
        ffrs.CIRC16(32, 16, 4, 2, simd_x16=False, simd_x8=False),
        ffrs.CIRC16(128, 8, 128, 128 // 16, simd_x16=False, simd_x8=False),
        ffrs.CIRC16(512, 4, 512, 512 // 16, simd_x16=False, simd_x8=False),
        # ffrs.CIRC(1024, 4, 1024, 1024 // 16, simd_x16=False, simd_x8=False),
        # ffrs.CIRC(4096, 4, 4096, 4096 // 16, simd_x16=False, simd_x8=False),
        ffrs.CIRC16(4, 2, 4, 2, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.CIRC16(8, 2, 8, 2, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.CIRC16(16, 4, 16, 4, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.CIRC16(32, 16, 4, 2, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.CIRC16(128, 8, 128, 128 // 16, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.CIRC16(512, 4, 512, 512 // 16, simd_x16=False, simd_x8=False, simd_x4=False),
        # ffrs.CIRC(1024, 4, 1024, 1024 // 16, simd_x16=False, simd_x8=False, simd_x4=False),
        # ffrs.CIRC(4096, 4, 4096, 4096 // 16, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.CIRC16(128, 2, 32, 2, 32),
        ffrs.CIRC16(512, 4, 128, 8, 8),
        ffrs.CIRC16(128, 8, 128, 128 // 16, 3),
        ffrs.CIRC16(128, 8, 128, 128 // 16, 4),
        ffrs.CIRC16(128, 8, 128, 128 // 16, 4, simd_x16=False),
        ffrs.CIRC16(128, 8, 128, 128 // 16, 4, simd_x16=False, simd_x8=False),
        ffrs.CIRC16(128, 8, 128, 128 // 16, 4, simd_x16=False, simd_x8=False, simd_x4=False),
    ],
)
class TestCircPower2(BaseTestCIRC):
    pass


@pytest.mark.parametrize(
    "rs",
    [
        ffrs.CIRC16(4 * 3, 2, 8 * 3, 2),
        ffrs.CIRC16(100, 2, 100, 2),
        ffrs.CIRC16(1026, 2, 128 * 3, 128),
        ffrs.CIRC16(100, 2, 100, 2, 3),
        ffrs.CIRC16(1026, 2, 128 * 3, 128, 7),
    ],
)
class TestCircMult(BaseTestCIRC):
    pass


@pytest.fixture(scope="class")
def setup_logger():
    level = ffrs.par.log.level
    ffrs.par.log.setLevel(logging.INFO)
    ffrs.set_logger(ffrs.par.log)
    yield
    ffrs.set_logger(None)
    ffrs.par.log.setLevel(level)


@pytest.mark.parametrize(
    "rs",
    [ffrs.CIRC16(256, 2, 512, 256, 128)],
)
@pytest.mark.usefixtures("setup_logger")
class TestCircSingle:
    pass
