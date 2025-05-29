import itertools
import pytest
import random

import ffrs

GF = ffrs.GFi32(65537, 3)

random.seed(42)

def randbytes(n, start=0, stop=256):
    return bytearray(random.randrange(start=start, stop=stop) for _ in range(n))

def field_samples(start=0):
    samples = 1000
    n = 16

    # Yield first n elements
    for i in range(start, min(n, GF.field_elements)):
        yield i

    # Yield last n elements
    for i in range(max(n, GF.field_elements - n), GF.field_elements):
        yield i

    for _ in range(samples - 2 * n):
        yield random.randrange(start, GF.field_elements)

@pytest.mark.parametrize('fn, id', [
    (GF.add, 0),
    (GF.sub, 0),
    (GF.mul, 1),
    (GF.div, 1),
])
def test_scalar_identity(fn, id):
    for a in field_samples():
        assert fn(a, id) == a, a


@pytest.mark.parametrize('fn1, fn2', [
    (GF.exp, GF.log),
])
def test_scalar_inverse1(fn1, fn2):
    for a in field_samples(start=1):
        assert fn2(fn1(a)) == a, a


@pytest.mark.parametrize('fn1, fn2', [
    (GF.add, GF.sub),
    (GF.mul, GF.div),
])
def test_scalar_inverse2(fn1, fn2):
    for a, b in itertools.product(field_samples(), field_samples(start=1)):
        res = fn1(a, b)
        res_inv = fn2(res, b)
        assert res_inv == a, (a, b)


@pytest.mark.parametrize('fn', [
    GF.add,
    GF.mul,
])
def test_scalar_commutativity(fn):
    for a, b in itertools.product(field_samples(), field_samples()):
        assert fn(a, b) == fn(b, a)


def test_scalar_inv_div():
    for a in field_samples(start=1):
        assert GF.inv(a) == GF.div(1, a) != 0


def test_scalar_inv_inv():
    for a in field_samples(start=1):
        assert GF.inv(a) == GF.div(1, a) != 0


def test_scalar_exp_pow():
    for a in field_samples():
        assert GF.exp(a) == GF.pow(GF.primitive, a) != 0

