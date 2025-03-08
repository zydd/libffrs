import pytest

from ffrs.reference import GF
from ffrs.reference.linalg import *

GF256 = GF(2, 8)


@pytest.mark.parametrize('n', [1, 2, 3, 4, 5, 10, 20])
def test_matmul(n):
    id1 = identity(GF256, n)
    id2 = identity(GF256, n)
    id3 = identity(GF256, n)

    assert id1 == matmul(GF256, id2, id3)


@pytest.mark.parametrize('n', [1, 2, 3, 4, 5, 10, 20])
def test_vandermonde_inverse(n):
    id = identity(GF256, n)
    v = vandermonde_ntt(GF256(2), n)
    vi = inverse(GF256, v)

    assert v == vandermonde_ntt(GF256(2), n)
    assert id == matmul(GF256, v, vi)
    assert inverse(GF256, vi) == v
