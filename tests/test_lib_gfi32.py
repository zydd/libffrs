import itertools
import pytest
import random

import ffrs
import ffrs.reference
import ffrs.reference.ntt as ref_ntt


GF = ffrs.GFi32(65537, 3)
GF_ref = lambda x: ffrs.reference.F(GF.prime, x)
random.seed(42)


def to_int_list(buf, size, byteorder="little"):
    return [int.from_bytes(buf[i:i+size], byteorder) for i in range(0, len(buf), size)]


def from_int_list(lst, size, byteorder="little"):
    buf = bytearray(len(lst) * size)
    for i, v in enumerate(lst):
        buf[i * size:(i + 1) * size] = int.to_bytes(v, size, byteorder)
    return buf


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


def test_ntt():
    size = 8
    w = next(i for i in range(2, GF.prime) if i**size % GF.prime == 1)
    buf = [random.randrange(0, GF.prime) for _ in range(size)]

    res_buf = GF.ntt(from_int_list(buf, 4), w)
    res = to_int_list(res_buf, 4)

    ref = [int(x) for x in ref_ntt.ntt(GF_ref, GF_ref(w), buf)]

    print()
    print("w:",w )
    print("input:      ", buf)
    print("reference:  ", ref)
    print("result:     ", res)

    assert ref == res
