import math
import random
import ffrs.reference as ref
import ffrs.reference.ntt as ref_ntt
import ffrs.reference.rs as ref_rs
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted

GF = ref.GF(65537, 1, 3)
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

    while l < max(deg, 4):
        print("inverse_modulo", l)

        a = a * (P([2]) - a * h)
        A = [a * (GF(2) - a * h) for a, h in zip(A, H)]

        # a = a % x**deg
        a.x = a.x[:deg]

        print(a)

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
    rev_g = P(g.x[::-1])

    mod = f.deg() - g.deg() + 1
    print("rev_g", rev_g.x)
    rev_g_i = inverse_modulo(rev_g, mod)
    print("rev_g_i", rev_g_i.x)
    rev_rev_g_i = P(rev_g_i.x[::-1])
    print("rev_rev_g_i", rev_rev_g_i.x)
    q = f * rev_rev_g_i
    print("f * rev_rev_g_i", q.x)
    q = P(q.x[-mod:])

    return q


def test():
    f = P([256, 32897, 0, 32896])
    g = P([0, 255, 256])
    q = f // g  # 65281, 32897
    r = f % g

    assert g * q == ntt_mul(g, q), (g * q, ntt_mul(g, q))

    print(f"f = gq+r = {f.x}")
    print(f"g = {g.x}")
    print(f"q = {q.x}")
    print(f"r = {r.x}")

    # coef = g.x[-1]
    # g = P([a // coef for a in g.x])
    # f = P([a // coef for a in f.x])
    rec_q = ntt_div(f, g)

    print(f"f/g = {rec_q.x}")
    # print(f"f%g = {(f - rec_q * g) * P([coef])}")
    print()

    assert q == rec_q, (q, rec_q)


if __name__ == "__main__":
    test()
