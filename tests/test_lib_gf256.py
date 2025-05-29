
# test_lib_gf256.py
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

import functools
import itertools
import ffrs.reference.ntt
import pytest
import random

import ffrs


GF256 = ffrs.GF256()
GF256_ref = ffrs.reference.GF(2, 8, GF256.primitive, GF256.poly1 | 0x100)

random.seed(42)

def randbytes(n, start=0, stop=256):
    return bytearray(random.randrange(start=start, stop=stop) for _ in range(n))


@pytest.mark.parametrize('fn, id', [
    (GF256.add, 0),
    (GF256.sub, 0),
    (GF256.mul, 1),
    (GF256.div, 1),
])
def test_scalar_identity(fn, id):
    for a in range(256):
        assert fn(a, id) == a, a


@pytest.mark.parametrize('fn1, fn2', [
    (GF256.exp, GF256.log),
])
def test_scalar_inverse1(fn1, fn2):
    for a in range(1, 256):
        assert fn2(fn1(a)) == a, a


@pytest.mark.parametrize('fn1, fn2', [
    (GF256.add, GF256.sub),
    (GF256.mul, GF256.div),
])
def test_scalar_inverse2(fn1, fn2):
    for a, b in itertools.product(range(256), range(1, 256)):
        assert fn2(fn1(a, b), b) == a, (a, b)


@pytest.mark.parametrize('fn', [
    GF256.add,
    GF256.sub,
    GF256.mul,
])
def test_scalar_commutativity(fn):
    for a, b in itertools.product(range(256), range(256)):
        assert fn(a, b) == fn(b, a)


def test_scalar_inv_div():
    for a in range(1, 256):
        assert GF256.inv(a) == GF256.div(1, a) != 0


def test_scalar_inv_inv():
    for a in range(1, 256):
        assert GF256.inv(a) == GF256.div(1, a) != 0


def test_scalar_exp_pow():
    for a in range(256):
        assert GF256.exp(a) == GF256.pow(GF256.primitive, a) != 0


# def test_pybind_u8_type_error():
#     # TODO: open issue on pybind11 due to cryptic error message
#     ffrs.GF256(poly1=0x11d)  # expects uint8_t, fails because 0x11d > 255


def test_poly_add_identity():
    for size in range(32):
        a = randbytes(size)
        b = bytearray(size)
        res = GF256.poly_add(a, b)

        assert res == a


@pytest.mark.parametrize('fn', [
    GF256.poly_add,
    GF256.poly_sub,
])
def test_poly_different_sizes(fn):
    for size_a, size_b in itertools.product(range(32), range(32)):
        a = randbytes(size_a)
        b = randbytes(size_b)
        size_max = max(size_a, size_b)
        res = fn(a, b)
        res_trunc = fn(a[:size_max], b[:size_max])

        assert len(res) == size_max

        if len(res) > 0:
            assert bytes(res).endswith(bytes(res_trunc))


def test_poly_add_sub():
    for size in range(32):
        a = randbytes(size)
        b = randbytes(size)
        assert GF256.poly_sub(GF256.poly_add(a, b), b) == a, (a, b)


def test_poly_divmod_mul_add():
    for size_a in range(32):
        for size_b in range(1, size_a):
            a = randbytes(size_a)
            b = randbytes(size_b)

            # FIXME: test fails if b[0] == 0
            b[0] |= 1

            res_q, res_r = GF256.poly_divmod(a, b)

            assert a == GF256.poly_add(GF256.poly_mul(res_q, b), res_r), (a, b)


def test_poly_divmod_mod():
    for size_a in range(32):
        for size_b in range(1, size_a):
            a = randbytes(size_a)
            b = randbytes(size_b)

            res_divmod_q, res_divmod_r = GF256.poly_divmod(a, b)
            res_mod_r = GF256.poly_mod(a, b)

            assert res_divmod_r == res_mod_r


def test_poly_mod_x_n():
    for size_a in range(32):
        for size_b in range(4, size_a):
            a = randbytes(size_a)
            b = randbytes(size_b)
            b[0] = 1

            x_n = bytearray(size_b)
            x_n[0] = 1

            res_r = GF256.poly_mod(GF256.poly_mul(a, x_n), b)

            res = GF256.poly_mod_x_n(a, b[1:])

            assert res_r == res

def test_poly_eval():
    for size in range(1, 32):
        roots = randbytes(size, start=1)

        # Construct a polynomial with roots r_i: ∏(1 - x/r_i)
        a = bytearray([1])
        for i in range(size):
            a = GF256.poly_mul(a, bytearray([GF256.inv(GF256.sub(0, roots[i])), 1]))

        for i in range(size):
            assert GF256.poly_eval(a, roots[i]) == 0


def test_mul8():
    for _ in range(256):
        a = randbytes(8)
        b = randbytes(8)

        ref = bytearray(map(lambda a: GF256.mul(*a), zip(a, b)))

        assert GF256.mul8(a, b) == ref


def test_poly_eval8():
    for _ in range(32):
        a = randbytes(8)
        b = bytearray(int(GF256_ref.gen(i)) for i in range(len(a)))

        res = list(GF256.poly_eval8(a, b))
        ref = [int(x) for x in ffrs.reference.ntt.ntt(GF256_ref, GF256_ref.a, a[::-1])]

        assert res == ref
