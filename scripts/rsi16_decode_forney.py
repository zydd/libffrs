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


x = ref.P(GF, [0, 1])


def sugiyama(synds):
    R2 = x ** len(synds)
    R1 = ref.P(GF, (synds))
    A2 = ref.P(GF, [0])
    A1 = ref.P(GF, [1])
    iteration = 0
    while R1.deg() >= len(synds) // 2:
        print(f"\nIteration: {iteration}")
        print("R2:", list(map(int, R2.x)))
        print("R1:", list(map(int, R1.x)))
        print("A2:", list(map(int, A2.x)))
        print("A1:", list(map(int, A1.x)))
        iteration += 1

        Q = R2 // R1

        print("Q:", list(map(int, Q.x)))
        breakpoint()
        print("---")
        print("R1 * Q:", list(map(int, (Q * R1).x)))
        print("A1 * Q:", list(map(int, (Q * A1).x)))
        print("---")

        t = A2 - Q * A1
        print("A2 - A1 * Q:", list(map(int, t.x)))
        A2 = A1
        A1 = t

        t = R2 - Q * R1
        print("R2 - R1 * Q:", list(map(int, t.x)))
        R2 = R1
        R1 = t

    print()
    print("R2:", list(map(int, R2.x)))
    print("R1:", list(map(int, R1.x)))
    print("A2:", list(map(int, A2.x)))
    print("A1:", list(map(int, A1.x)))
    print()

    locator = ref.P(GF, [a // GF(A1.x[0]) for a in A1.x])
    evaluator = ref.P(GF, [a // GF(A1.x[0]) for a in R1.x])
    print("locator:", list(map(int, locator.x)))
    print("evaluator:", list(map(int, evaluator.x)))
    return locator, evaluator


def forney(locator, evaluator):
    err_pos = [i for i in range(size) if locator.eval(w**-i) == 0]

    err_val = []
    for pos in err_pos:
        err_val.append(-(w**pos) * evaluator.eval(w**-pos) // locator.deriv().eval(w**-pos))

    err_pos = [rbo(size, i) for i in err_pos]
    print("err_pos:", err_pos)
    print("err_val:", to_int_list(err_val))
    return err_pos, err_val


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
    err_pos = [random.randint(0, size - 1) for _ in range(1)]
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

    locator, evaluator = sugiyama(E[:ecc])

    if L != locator:
        print("**** sugiyama failed ****")

    forney_err_pos, forney_err_val = forney(locator, evaluator)
    assert set(err_pos) == set(forney_err_pos)
    assert all(e[pos] == forney_err_val[i] for i, pos in enumerate(forney_err_pos))

    print()
    print(64 * "-")
    synd = [256, 32897, 0, 32896]
    # synd = [128, 128, 128, 128, ]
    synd = list(map(GF, synd))
    sugiyama(synd)


if __name__ == "__main__":
    test()
