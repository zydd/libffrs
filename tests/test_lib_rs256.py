
#  test_lib_rs.py
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

import pytest
import random

import ffrs


random.seed(42)

def randbytes(n):
    return bytearray(random.randrange(256) for _ in range(n))


@pytest.mark.parametrize('rs', [
    ffrs.RS256(ecc_len=l) for l in
        list(range(1, 16 + 1)) + [31, 32, 33, 63, 64, 65, 127, 128, 129, 253, 254]])
class TestRS:
    def test_encode(self, rs):
        for size in range(rs.gf.field_elements-1 - rs.ecc_len):
            a = randbytes(size)
            rem = rs.gf.poly_mod_x_n(a, rs.generator[1:])

            msg_a = a + bytearray(rs.ecc_len)
            rs.encode(msg_a)

            assert msg_a == a + rem


    def test_encode_decode(self, rs):
        for size in range(rs.gf.field_elements-1 - rs.ecc_len):
            a = randbytes(size)

            msg_a = a + bytearray(rs.ecc_len)
            rs.encode(msg_a)

            msg_b = bytearray(msg_a)

            rs.decode(msg_b)
            assert msg_a == msg_b

            error_positions = set()
            while len(error_positions) < rs.ecc_len//2:
                i = random.randrange(len(msg_b))
                if i in error_positions:
                    continue
                error_positions.add(i)
                msg_b[i] ^= random.randrange(1, 256)

            assert rs.decode(msg_b) is True
            assert msg_a == msg_b


    def test_decode_fail(self, rs):
        for size in range(rs.gf.field_elements-1 - rs.ecc_len):
            a = randbytes(size)

            msg_a = a + bytearray(rs.ecc_len)
            rs.encode(msg_a)

            msg_b = bytearray(msg_a)

            rs.decode(msg_b)
            assert msg_a == msg_b

            error_positions = set()
            while len(error_positions) < rs.ecc_len//2 + 1:
                i = random.randrange(len(msg_b))
                if i in error_positions:
                    continue
                error_positions.add(i)
                msg_b[i] ^= random.randrange(1, 256)

            dec_status = rs.decode(msg_b)

            if dec_status is True:
                # True status here means the corrupted message has the same polynomial roots as the original
                msg_c = msg_b[:-rs.ecc_len] + bytearray(rs.ecc_len)
                rs.encode(msg_c)
                assert msg_b == msg_c

            assert msg_a != msg_b


    def test_encode_blocks(self, rs):
        for size in range(1, rs.gf.field_elements-1 - rs.ecc_len):
            a = randbytes(size)
            msg_a = a + bytearray(rs.ecc_len)
            rs.encode(msg_a)

            # block size == input size
            assert msg_a == rs.encode_blocks(a, len(a))

            # blocks size > input size
            max_block_size = 254 - rs.ecc_len
            assert msg_a == rs.encode_blocks(a, min(len(a) + 1, max_block_size))
            assert msg_a == rs.encode_blocks(a, min(len(a) + 10, max_block_size))


    def test_encode_blocks_empty(self, rs):
        assert bytearray() == rs.encode_blocks(bytearray(), 0)
        assert bytearray() == rs.encode_blocks(bytearray(10), 0)
        assert bytearray() == rs.encode_blocks(bytearray(), 10)


    def test_encode_blocks_lenghts(self, rs):
        for data_size in [1, 2, 3, 200, 201, 254]:
            for block_size in range(1, data_size + 1):
                a = randbytes(data_size)

                msg = rs.encode_blocks(a, block_size)

                for i in range(data_size // block_size):
                    msg_i = a[i * block_size:][:block_size] + bytearray(rs.ecc_len)
                    rs.encode(msg_i)

                    assert msg_i == msg[i * (block_size + rs.ecc_len):][:block_size + rs.ecc_len]

                tail_size = data_size % block_size
                if tail_size > 0:
                    msg_i = a[-tail_size:] + bytearray(rs.ecc_len)
                    rs.encode(msg_i)
