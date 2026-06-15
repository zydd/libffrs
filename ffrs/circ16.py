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

from .util import b64hex, b64hex_dec


class CIRC16(libffrs.CIRC16):
    __doc__ = libffrs.CIRC16.__doc__

    def __repr__(self):
        return (
            f"CIRC("
            f"{self.inner_block_len}, {self.inner_ecc_len}"
            f", {self.outer_block_len}, {self.outer_ecc_len}"
            f", {self.interleave})"
        )

    def serialize(self):
        """
        Serialize CIRC16 config

        Example:
        circ16 ip:g01 ig:3 ib:x4 im:x0 ie:4 op:g01 og:3 ob:80 om:40 oe:40 i:1
        """
        return (
            "circ16"
            f" ip:{b64hex(self.rsi.gf.field_elements)} ig:{b64hex(self.rsi.gf.primitive)} ib:{b64hex(self.inner_block_len)} im:{b64hex(self.inner_message_len)} ie:{b64hex(self.inner_ecc_len)}"
            f" op:{b64hex(self.rso.gf.field_elements)} og:{b64hex(self.rso.gf.primitive)} ob:{b64hex(self.outer_block_len)} om:{b64hex(self.outer_message_len)} oe:{b64hex(self.outer_ecc_len)}"
            f" i:{b64hex(self.interleave)}"
        ).encode("ascii")

    @staticmethod
    def deserialize(buf):
        """
        Create CIRC16 instance from serialized form
        """
        string = buf.decode("ascii")
        assert string.startswith("circ16 ")
        params = {}
        for part in string[7:].split():
            key, value = part.split(":", 1)
            params[key] = b64hex_dec(value)

        assert params["ip"] == params["op"] == 0x10001
        assert params["ig"] == params["og"]
        assert params["ib"] == params["im"] + params["ie"]
        assert params["ob"] == params["om"] + params["oe"]

        return CIRC16(
            inner_block_len=params["ib"],
            inner_ecc_len=params["ie"],
            outer_block_len=params["ob"],
            outer_ecc_len=params["oe"],
            interleave=params["i"],
            primitive=params["ig"],
        )
