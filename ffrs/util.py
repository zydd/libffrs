#  util.py
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


_b64hex_alphabet = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+-"


def b64hex(val: int) -> str:
    result = ""
    while val:
        result = _b64hex_alphabet[val & 0x3F] + result
        val >>= 6
    return result


def b64hex_dec(val: str) -> int:
    result = 0
    for c in val:
        result <<= 6
        result |= _b64hex_alphabet.index(c)
    return result
