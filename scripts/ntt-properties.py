import random

import ffrs.reference as ref
import ffrs.reference.ntt as ref_ntt
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted

GF = ref.GF(257, 1, 5)
P = lambda a: ref.P(GF, a)

size = 8

w = GF.exp(GF.log(1) // GF(size))
print(f"w: {w}")
print()
assert w**size == 1

ntt = lambda arr: ref_ntt.ntt(GF, w, arr)
intt = lambda arr: ref_ntt.intt(GF, w, arr)

fmt = lambda arr: [a if a <= 257 // 2 else a - 257 for a in (int(a) for a in arr)]
rpad = lambda arr: arr + [GF(0)] * (size - len(arr))


def time_shift(Y, y_shift, shift):
    return [a * w ** (shift * i) for i, a in enumerate(Y)], y_shift + shift


def poly_add(P, p_shift, Q, q_shift):
    Qs, _qs_shift = time_shift(Q, q_shift, p_shift - q_shift)
    return [p + q for p, q in zip(P, Qs)], p_shift


def poly_mul(P, p_shift, Q, q_shift):
    return [p * q for p, q in zip(P, Q)], p_shift + q_shift - 1


def test_properties():
    y = [GF(random.randint(0, 256)) for i in range(size)]
    Y = ntt(y)

    print(f"y: {fmt(y)}")
    print(f"Y: {fmt(Y)}")
    print()
    assert intt(Y) == y, intt(Y)

    # Shift
    shift = 3 % size
    y_shift = y[shift:] + y[:shift]
    Y_shift = [a * w ** (-shift * i) for i, a in enumerate(Y)]

    print(f"y_shift: {fmt(y_shift)}")
    print(f"Y_shift: {fmt(Y_shift)}")
    print()
    assert ntt(y_shift) == Y_shift

    # Reverse
    y_rev = y[::-1]
    Y_rev = [a * w**-i for i, a in enumerate(Y[:1] + Y[:0:-1])]

    print(f"y_rev: {fmt(y_rev)}")
    print(f"Y_rev: {fmt(Y_rev)}")
    print(f"Y_rev: {fmt(ntt(y_rev))}")
    print()
    assert ntt(y_rev) == Y_rev

    # print(f"ratio: {fmt([a // b if int(b) != 0 else -666 for a, b in zip(Y_rev, ntt(y_rev))])}")
    # print(f"ratio logw: {fmt([GF.log(w, a // b) if int(b) != 0 else -666 for a, b in zip(Y_rev, ntt(y_rev))])}")


def test_add_mul():
    a = [1, 2, 3]
    b = [4, 5]
    c = [6, 7]

    a, b, c = P(a), P(b), P(c)

    print(f"a:      {a}")
    print(f"b:      {b}")
    print(f"c:      {c}")
    print(f"a * b:  {a * b}")
    print(f"a*b+c:  {a * b + c}")
    print()

    A = ntt(rpad(a.x[::-1]))
    B = ntt(rpad(b.x[::-1]))
    C = ntt(rpad(c.x[::-1]))
    a_shift = len(a.x)
    b_shift = len(b.x)
    c_shift = len(c.x)

    print(f"A:      {fmt(A)} ({a_shift})")
    print(f"B:      {fmt(B)} ({b_shift})")
    print(f"B:      {fmt(C)} ({c_shift})")

    F = [a * b for a, b in zip(A, B)]
    f_shift = a_shift + b_shift - 1
    print(f"A(*)B:  {fmt(F)} ({f_shift})")
    print(f"a * b:  {fmt(intt(F))} ({f_shift})")
    print()

    S, s_shift = time_shift(C, c_shift, f_shift - c_shift)
    print(f"C (sh): {fmt(S)} ({s_shift})")
    print(f"c (sh): {fmt(intt(S))} ({s_shift})")
    print()

    assert f_shift == s_shift

    R = [f + s for f, s in zip(F, S)]
    r_shift = f_shift
    print(f"A(*)B+C:{fmt(R)} ({s_shift})")
    print(f"a*b+c:  {fmt(intt(R))} ({s_shift})")
    print()

    assert rpad((a * b + c).x[::-1]) == intt(R)

    R2, r2_shift = poly_add(*poly_mul(A, a_shift, B, b_shift), C, c_shift)
    print(f"A(*)B+C:{fmt(R2)} ({r2_shift})")
    print(f"a*b+c:  {fmt(intt(R2))} ({r2_shift})")

    assert R2 == R
    assert r2_shift == r_shift


if __name__ == "__main__":
    # test_properties()
    test_add_mul()
