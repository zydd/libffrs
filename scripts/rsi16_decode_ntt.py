import math
import random
import ffrs.reference as ref
import ffrs.reference.ntt as ref_ntt
import ffrs.reference.rs as ref_rs
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted


GF = ref.GF(0x10001, 1, 3)
ntt = lambda w, x: ref_ntt.ntt(GF, w, x)
intt = lambda w, x: ref_ntt.intt(GF, w, x)


roots_of_unity = [1, 1, 256, 4096, 64, 8, 8224, 13987, 282, 15028, 19139, 3668, 11088, 6561, 81, 9]
# roots_of_unity ** 2 = [1, 1, 1, 256, 4096, 64, 8, 8224, 13987, 282, 15028, 19139, 3668, 11088, 6561, 81]


size = 16
ecc = 4
ecc_decimation = size // ecc
w = GF(roots_of_unity[round(math.log2(size))])
assert w**size == 1


def test():
    m = [GF(0)] * ecc + [GF(i + 1000) for i in range(size - ecc)]
    print("m:", to_int_list(m))

    M = ntt(w, rbo_sorted(m))

    print("M:", to_int_list(M))
    # print("m:", to_int_list(rbo_sorted(intt(w, M))))
    # M2 = ref_ntt.matmul(GF, ref_ntt.vandermonde_ntt(w, size), rbo_sorted(m))
    # print("V:", to_int_list(M2))
    print()

    print("C = M[:ecc] * (size // ecc)")
    C = M[:ecc] * ecc_decimation
    print("C:", to_int_list(C))

    c = rbo_sorted(intt(w, C))
    print("c:", to_int_list(c))

    print()

    print("t = m - c")
    t = [m[i] - c[i] for i in range(size)]
    print("t:", to_int_list(t))
    T = ntt(w, rbo_sorted(t))
    print("T:", to_int_list(T))
    print()

    e = [GF(0) for i in range(size)]
    err_pos = [random.randint(0, size-1) for _ in range(ecc)]
    for pos in err_pos:
        e[pos] += GF(random.randint(0, 0x10000))
    print("e:", to_int_list(e))
    E = ntt(w, rbo_sorted(e))
    E[ecc:] = [0] * (size - ecc)
    print("E:", to_int_list(E))
    print()

    print("r = t + e")
    r = [t[i] + e[i] for i in range(size)]
    print("r:", to_int_list(r))

    R = ntt(w, rbo_sorted(r))
    print("R:", to_int_list(R))
    print()

    L = ref.P(GF, [1])
    L = ref_rs.locator(GF, w, [rbo(size, pos) for pos in err_pos])

    print("L:", L)
    # print("L(err_pos):", [int(L.eval(w**-rbo(size, pos))) for pos in err_pos])
    print()

    for j in range(ecc, size):
        E[j] = -sum((l * e for l, e in zip(L.x[:0:-1], E[j - len(err_pos):j])), start=GF(0))

    M2 = [R[j] - E[j] for j in range(size)]
    m2 = rbo_sorted(intt(w, M2))
    print("M2:", to_int_list(M2))
    print("m2:", to_int_list(m2))
    print()
    assert t == m2


    e2 = rbo_sorted(intt(w, E))

    print("E:", to_int_list(E))
    print("e2:", to_int_list(e2))

    assert e == e2


if __name__ == "__main__":
    test()
