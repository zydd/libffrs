
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

from ffrs.rsi16md import RSi16md


class CIRC:
    def __init__(self, inner_block_len, inner_ecc_len, outer_block_len, outer_ecc_len):
        self.rsi = RSi16md(inner_block_len, ecc_len=inner_ecc_len)
        self.rso = RSi16md(outer_block_len, ecc_len=outer_ecc_len, interleave=self.rsi.message_len)

    def encode(self, data):
        assert len(data) % self.message_len == 0
        inner = self.rsi.encode_blocks(data)
        outer = self.rso.encode_chunk(data)
        return inner + outer

    @property
    def message_len(self):
        return self.rsi.message_len * self.rso.message_len

    @property
    def ecc_len(self):
        return self.rsi.ecc_len * self.rso.message_len + self.rso.ecc_len * self.rso.interleave

    def check(self, data, ecc, synd=None):
        synd = synd or bytearray(self.ecc_len)

        assert len(data) % self.message_len == 0
        assert len(ecc) % self.ecc_len == 0
        assert len(data) // self.message_len == len(ecc) // self.ecc_len
        assert len(synd) == len(ecc)

        inner_ecc_len = self.rsi.ecc_len * self.rso.message_len

        inner_ok = self.rsi.check_blocks(data, ecc[:inner_ecc_len], synd[:inner_ecc_len])
        outer_ok = self.rso.check_chunk(data, ecc[inner_ecc_len:], synd[inner_ecc_len:])
        return inner_ok and outer_ok, synd
    
    def repair(self, data, ecc, synd=None):
        if not synd:
            synd = bytearray(self.ecc_len)
            integrity_ok, synd = self.check(data, ecc, synd)
            if integrity_ok:
                return data

        inner_ecc_len = self.rsi.ecc_len * self.rso.message_len

        nok_rows = self._nok_rows(synd)
        if len(nok_rows) <= self.rso.ecc_len:
            self.rso.repair_chunk_error_locations(data, nok_rows)
            return data
        else:
            self.rsi.repair_blocks(data, synd[:inner_ecc_len])
            # TODO: find corrupted block start if len(nok_rows) == self.rso.ecc_len + 1
            # TODO: retry outer repair after

        raise NotImplementedError
