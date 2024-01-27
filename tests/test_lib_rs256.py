
#  test_lib_rs256.py
#
#  Copyright 2024 Gabriel Machado
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

import random

import ffrs

RS256 = ffrs.RS256(ecc_len=8)


random.seed(42)

def randbytes(n):
    return bytearray(random.randrange(256) for _ in range(n))


def test_encode():
    for size in range(255 - RS256.ecc_len):
        a = randbytes(size)
        rem = RS256.gf.poly_mod_x_n(a, RS256.generator[1:])

        msg_a = a + bytearray(RS256.ecc_len)
        RS256.encode(msg_a)

        assert msg_a == a + rem


def test_encode_decode():
    for size in range(255 - RS256.ecc_len):
        a = randbytes(size)
        rem = RS256.gf.poly_mod_x_n(a, RS256.generator[1:])

        msg_a = a + bytearray(RS256.ecc_len)
        RS256.encode(msg_a)

        msg_b = bytearray(msg_a)

        RS256.decode(msg_b)
        assert msg_a == msg_b

        error_positions = set()
        while len(error_positions) < RS256.ecc_len//2:
            i = random.randrange(len(msg_b))
            if i in error_positions:
                continue
            error_positions.add(i)
            msg_b[i] ^= random.randrange(1, 256)

        assert msg_a != msg_b
        assert RS256.decode(msg_b) is True
        assert msg_a == msg_b


def test_decode_fail():
    for size in range(255 - RS256.ecc_len):
        a = randbytes(size)

        msg_a = a + bytearray(RS256.ecc_len)
        RS256.encode(msg_a)

        msg_b = bytearray(msg_a)

        RS256.decode(msg_b)
        assert msg_a == msg_b

        error_positions = set()
        while len(error_positions) < RS256.ecc_len//2 + 1:
            i = random.randrange(len(msg_b))
            if i in error_positions:
                continue
            error_positions.add(i)
            msg_b[i] ^= random.randrange(1, 256)

        dec_status = RS256.decode(msg_b)

        if dec_status is True:
            # True status here means the corrupted message has the same polynomial roots as the original
            msg_c = msg_b[:-RS256.ecc_len] + bytearray(RS256.ecc_len)
            RS256.encode(msg_c)
            assert msg_b == msg_c

        assert msg_a != msg_b


def test_encode_blocks():
    for data_size in [1, 2, 3, 200, 201, 254]:
        for block_size in range(1, data_size + 1):
            a = randbytes(data_size)

            msg = RS256.encode_blocks(a, block_size)

            for i in range(data_size // block_size):
                msg_i = a[i * block_size:][:block_size] + bytearray(RS256.ecc_len)
                RS256.encode(msg_i)

                assert msg_i == msg[i * (block_size + RS256.ecc_len):][:block_size + RS256.ecc_len]

            tail_size = data_size % block_size
            if tail_size > 0:
                msg_i = a[-tail_size:] + bytearray(RS256.ecc_len)
                RS256.encode(msg_i)
