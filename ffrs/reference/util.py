
#  util.py
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
import random

import ffrs.reference


def sample_field(GF, start=0, samples=1000):
    if GF.field_elements < samples:
        return range(GF.field_elements)

    n = 8

    # Yield first n elements
    for i in range(start, min(n, GF.field_elements)):
        yield i

    # Yield last n elements
    for i in range(max(n, GF.field_elements - n), GF.field_elements):
        yield i

    for _ in range(samples - 2 * n):
        yield random.randrange(start, GF.field_elements)


def to_bytearray(v, sizeof, byteorder="little"):
    if type(v) in [bytearray, bytes, str]:
        return bytearray(v, "utf-8")

    buf = bytearray(len(v) * sizeof)
    for i, v in enumerate(v):
        buf[i * sizeof:(i + 1) * sizeof] = int.to_bytes(v, sizeof, byteorder)
    return buf


def to_int_list(v, sizeof=None, byteorder="little"):
    if type(v) is list:
        if len(v) == 0: return []

        if type(v[0]) in [ffrs.reference.F, ffrs.reference.Fp]:
            return list(map(int, v))
        else:
            raise TypeError
    else:
        assert sizeof, "Must specify sizeof for byte objects"
        return [int.from_bytes(v[i:i+sizeof], byteorder) for i in range(0, len(v), sizeof)]


def randbytes(n, start=0, stop=256):
    return bytearray(random.randrange(start=start, stop=stop) for _ in range(n))


def rbo(max, i):
    nbits = round(math.log2(max))
    return int(bin(i)[2:].zfill(nbits)[::-1], 2)


def rbo_sorted(a):
    assert len(a) and len(a) & (len(a) - 1) == 0, "len(a) must be a power of 2"
    return [a[rbo(len(a), i)] for i in range(len(a))]
