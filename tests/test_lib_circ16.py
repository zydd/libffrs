
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
    ffrs.CIRC(4096, 4, 4096, 4096 // 16),
])
class TestCIRC:
    def test_circ_encode(self, rs):
        buf = randbytes(rs.message_len)
        res = rs.encode(buf)
