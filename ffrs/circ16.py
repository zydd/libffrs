#  circ16.py
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

import libffrs


class CIRC(libffrs.CIRCi16):
    def check(self, data, ecc, synd=None):
        synd = synd or bytearray(self.ecc_size)

        assert len(data) % self.message_size == 0
        assert len(ecc) % self.ecc_size == 0
        assert len(data) // self.message_size == len(ecc) // self.ecc_size
        assert len(synd) == len(ecc)

        inner_ecc_size = self.rsi.ecc_size * self.rso.message_size

        inner_ok = self.rsi.check_blocks(data, ecc[:inner_ecc_size], synd[:inner_ecc_size])
        outer_ok = self.rso.check_chunk(data, ecc[inner_ecc_size:], synd[inner_ecc_size:])
        return inner_ok and outer_ok, synd

    def __repr__(self):
        return (
            f"CIRC("
            f"{self.inner_block_len}, {self.inner_ecc_len}"
            f", {self.outer_block_len}, {self.outer_ecc_len}"
            f", {self.outer_interleave})"
        )
