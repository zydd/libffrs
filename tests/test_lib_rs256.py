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

import itertools
import pytest
import random

import ffrs


def randbytes(n):
    return bytearray(random.randrange(256) for _ in range(n))


def test__init__():
    rs1 = ffrs.RS256(254, 222)
    rs2 = ffrs.RS256(ecc_len=32)
    rs3 = ffrs.RS256(254, 222, ecc_len=32)

    assert rs1.ecc_len == rs2.ecc_len == rs3.ecc_len == 32
    assert rs1.block_len == rs3.block_len == 254

    assert rs2.block_len == 255
    assert rs2.message_len == 255 - 32
    rs2.block_len = 48
    assert rs2.message_len == 48 - 32

    assert 222 == ffrs.RS256(254, 222).message_len
    assert 222 == ffrs.RS256(block_len=254, ecc_len=32).message_len
    assert 222 == ffrs.RS256(254, 222, ecc_len=32).message_len

    assert ffrs.RS256(message_len=8, ecc_len=9).block_len == 17
    assert ffrs.RS256(message_len=8, ecc_len=9).message_len == 8

    with pytest.raises(ValueError):
        ffrs.RS256(message_len=223)

    with pytest.raises(ValueError):
        ffrs.RS256(255, 223, ecc_len=33)  # 255 - 223 = 32

    with pytest.raises(ValueError):
        ffrs.RS256(254, 255)

    with pytest.raises(ValueError):
        ffrs.RS256(block_len=8, ecc_len=8)

    with pytest.raises(ValueError):
        ffrs.RS256(block_len=8, ecc_len=9)

    with pytest.raises(ValueError):
        ffrs.RS256()

    with pytest.raises(ValueError):
        ffrs.RS256(ecc_len=32).block_len = 32

    with pytest.raises(ValueError):
        ffrs.RS256(ecc_len=32).block_len = 0


@pytest.mark.parametrize(
    "rs",
    [
        ffrs.RS256(size, ecc_len=ecc)
        for size, ecc in itertools.product(
            range(255 - 16, 256, 16), list(range(1, 16 + 1)) + [31, 32, 33, 63, 64, 65, 127, 128, 129, 253, 254]
        )
        if size > ecc
    ],
)
class TestRS:
    def test_encode(self, rs):
        msg_a = randbytes(rs.message_size)
        rem = rs.gf.poly_mod_x_n(msg_a, rs.generator[1:])

        ecc = rs.encode(msg_a)

        assert ecc == rem

    def test_encode_decode(self, rs):
        msg_a = randbytes(rs.message_size)

        ecc_a = rs.encode(msg_a)

        msg_b = bytearray(msg_a)
        ecc_b = bytearray(ecc_a)

        rs.repair(msg_b, ecc_b)
        assert msg_a == msg_b
        assert ecc_a == ecc_b

        self._add_errors(msg_b, ecc_b, rs.ecc_len // 2)

        if rs.ecc_len >= 2:
            assert msg_a + ecc_a != msg_b + ecc_b

        assert rs.repair(msg_b, ecc_b) is True
        assert msg_a == msg_b
        assert ecc_a == ecc_b

    def test_decode_fail(self, rs):
        msg_a = randbytes(rs.message_size)
        ecc_a = rs.encode(msg_a)

        msg_b = bytearray(msg_a)
        ecc_b = bytearray(ecc_a)

        rs.repair(msg_b, ecc_b)
        assert msg_a == msg_b
        assert ecc_a == ecc_b

        self._add_errors(msg_b, ecc_b, rs.ecc_len // 2 + 1)

        dec_status = rs.repair(msg_b, ecc_b)

        if dec_status is True:
            # True status here means the corrupted message has the same polynomial roots as the original
            ecc_c = rs.encode(msg_b)
            assert ecc_b == ecc_c

        assert msg_a + ecc_a != msg_b + ecc_b

    def test__synds(self, rs):
        msg_a = randbytes(rs.message_size)
        ecc_a = rs.encode(msg_a)

        msg_b = bytearray(msg_a)
        ecc_b = bytearray(ecc_a)

        self._add_errors(msg_b, ecc_b, rs.ecc_len // 2)

        eval_gen_roots_a = bytearray(rs.gf.poly_eval(msg_a + ecc_a, x) for x in rs.generator_roots)
        eval_gen_roots_b = bytearray(rs.gf.poly_eval(msg_b + ecc_b, x) for x in rs.generator_roots)

        synd_a = rs._synds(msg_a + ecc_a)
        synd_b = rs._synds(msg_b + ecc_b)

        assert synd_a == bytearray(rs.ecc_len)
        assert synd_a == eval_gen_roots_a
        assert synd_b == eval_gen_roots_b

        if rs.ecc_len >= 2:
            assert msg_a + ecc_a != msg_b + ecc_b
            assert synd_b != synd_a

        rem_a = rs.gf.poly_mod(msg_a + ecc_a, rs.generator)
        rem_b = rs.gf.poly_mod(msg_b + ecc_b, rs.generator)

        assert synd_a == rem_a
        # assert synd_b == rem_b  # FIXME

    def test_encode_blocks_empty(self, rs):
        assert bytearray() == rs.encode(bytearray())

    def test_encode_blocks_multi(self, rs):
        for blocks in range(0, 5):
            for remainder in range(0, rs.message_size):
                data_size = blocks * rs.message_size + remainder
                msg = randbytes(data_size)

                ecc = rs.encode(msg)

                assert len(ecc) == rs.ecc_size * (blocks + int(remainder > 0))

                for i in range(data_size // rs.message_size):
                    msg_i = msg[i * rs.message_size:][:rs.message_size]
                    ecc_i = rs.encode(msg_i)

                    assert ecc_i == ecc[i * rs.ecc_size:][:rs.ecc_size]

                tail_size = data_size % rs.message_size
                if tail_size > 0:
                    msg_i = msg[-tail_size:]
                    ecc_i = rs.encode(msg_i)
                    assert ecc_i == ecc[-rs.ecc_size :]

    def _add_errors(self, msg, ecc, count):
        error_positions = set()
        while len(error_positions) < count:
            i = random.randrange(len(msg) + len(ecc))
            if i in error_positions:
                continue
            error_positions.add(i)

            if i < len(msg):
                msg[i] ^= random.randrange(1, 256)
            else:
                ecc[i - len(msg)] ^= random.randrange(1, 256)