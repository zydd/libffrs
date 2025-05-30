import pytest
import random

from ffrs.reference import *

from common import sample_field

fields = [
    GF(p) for p in [2, 3, 5, 7, 11, 13, 127, 257]] + [

    GF(2, 8, 2, 0x11d),
    GF(2, 1, 1),
    GF(3, 4),
    GF(5, 3),
    GF(7, 2),
    GF(11, 2),
    GF(2, 16, 2, 0x1002d)
]

@pytest.mark.parametrize("GF", fields)
def test_field(GF):
    zero = GF(0)
    one = GF(1)

    def assert_eq(a, b, l, r):
        assert l == r, (int(a), int(b), int(l), int(r))
        assert int(l) < GF.field_elements
        assert int(r) < GF.field_elements

    for i in sample_field(GF, start=1, samples=50):
        a = GF(i)
        assert i == int(a)

        try:
            a = int(a // GF(0))
            assert 0, "division by zero allowed"
        except ZeroDivisionError:
            pass

        assert_eq(a,    zero,   a * zero,       zero)
        assert_eq(zero, a,      zero // a,      zero)
        assert_eq(a,    one,    a * one,        a)
        assert_eq(a,    a,      a.inv() * a,    one)
        assert_eq(a,    a,      a.inv().inv(),  a)
        assert_eq(one,  a,      a.inv(),        one // a)

        for j in sample_field(GF, start=1, samples=50):
            b = GF(j)

            if hasattr(a, "_mul"):
                assert_eq(a, b, (a * b), (a._mul(b)))

            assert_eq(a, b, (a * b), (b * a))
            assert_eq(a, b, ((a * b) // b), (a))
            assert_eq(a, b, ((a * b) // a), (b))
            assert_eq(a, b, (a.inv() * b), (b // a))

            assert_eq(a, b, (a + b), (b + a))
            assert_eq(a, b, ((a + b) - a), (b))
            assert_eq(a, b, ((a + b) - b), (a))
            assert_eq(a, b, ((b - a) + a), (b))
            assert_eq(a, b, (-a + b), (b - a))


def test_gf256():
    GF256 = GF(2, 8)

    g = P(GF256, [1])
    i = GF256(1)
    for n in range(4):
        assert i == GF256.gen(n)

        g *= P(GF256, [i, 1])
        i *= GF256(2)
    assert str(g) == "𝑥⁴ + 15𝑥³ + 54𝑥² + 120𝑥 + 64"

    for i in range(4):
        assert int(g.eval(GF256.gen(i))) == 0, (i)

    assert len(list(GF.irr_polynomials(2, 8, 2))) == 16

