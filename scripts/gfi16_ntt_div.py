import math
import random
import ffrs.reference as ref
import ffrs.reference.ntt as ref_ntt
import ffrs.reference.rs as ref_rs
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted

GF = ref.GF(257, 1, 3)
P = lambda x: ref.P(GF, x)


ntt_poly_size = 16
for i in range(2, GF.p):
    ntt_poly_root = GF(i)
    if ntt_poly_root**ntt_poly_size == 1:
        break

ntt_poly = lambda x: ref_ntt.ntt(GF, ntt_poly_root, x.x + [GF(0)] * (ntt_poly_size - len(x.x)))
intt_poly = lambda x: (ref_ntt.intt(GF, ntt_poly_root, x))

x = P([0, 1])


def inverse_modulo(h, deg):
    # assert h % x == P([1])
    l = 1
    a = P([h.x[0].inv()])

    H = ntt_poly(h)
    A = ntt_poly(a)

    while l < deg:
        a = a * (P([2]) - a * h)
        A = [a * (GF(2) - a * h) for a, h in zip(A, H)]

        # a = a % x**deg
        a.x = a.x[:deg]

        l *= 2
        a = a
        A = A

    a2 = P(intt_poly(A)[:deg])

    assert a2 == a, (a, a2)

    return a


def ntt_mul(g, q):
    deg_f = g.deg() + q.deg() + 1

    rev_g = P(g.x[::-1])
    rev_q = P(q.x[::-1])

    rev_G = ntt_poly(rev_g)
    rev_Q = ntt_poly(rev_q)
    rev_F = [g * q for g, q in zip(rev_G, rev_Q)]

    G = ntt_poly(g)
    Q = ntt_poly(q)
    F = [g * q for g, q in zip(G, Q)]

    # print("rev_f", intt_poly(rev_F))
    # print("f    ", intt_poly(F))

    f1 = P(intt_poly(rev_F)[:deg_f][::-1])
    f2 = P(intt_poly(F)[:deg_f])
    assert f1 == f2
    return f1


def ntt_div(f, g):
    rev_f = P(f.x[::-1])
    rev_g = P(g.x[::-1])

    mod = x ** (f.deg() - g.deg() + 1)
    rev_g_i = inverse_modulo(rev_g, mod.deg())
    rev_q = rev_f * rev_g_i % mod

    q = P(rev_q.x[::-1]) * x ** (mod.deg() - rev_q.deg() - 1)
    return q


def test():
    g = P([13, 22, 66, 9, 13])
    q = P([0, 44, 28])
    r = P([11])
    f = g * q + r

    assert g * q == ntt_mul(g, q), (g * q, ntt_mul(g, q))

    print(f"g = {g}")
    print(f"q = {q}")
    print(f"r = {r}")
    print(f"f = gq+r = {f}")

    # coef = g.x[-1]
    # g = P([a // coef for a in g.x])
    # f = P([a // coef for a in f.x])
    rec_q = ntt_div(f, g)

    print(f"f/g = {rec_q}")
    # print(f"f%g = {(f - rec_q * g) * P([coef])}")
    print()

    assert rec_q == q, (q, rec_q)


if __name__ == "__main__":
    test()
