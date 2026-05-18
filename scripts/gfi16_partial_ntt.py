import math
import random
import ffrs.reference as ref
import ffrs.reference.ntt as ref_ntt
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted

GF = ref.GF(0x10001, 1, 3)

ntt = lambda w, x: ref_ntt.ntt(GF, w, rbo_sorted(x))
intt = lambda w, x: rbo_sorted(ref_ntt.intt(GF, w, x))


roots_of_unity = [1, 1, 256, 4096, 64, 8, 8224, 13987, 282, 15028, 19139, 3668, 11088, 6561, 81, 9]
# roots_of_unity ** 2 = [1, 1, 1, 256, 4096, 64, 8, 8224, 13987, 282, 15028, 19139, 3668, 11088, 6561, 81]


size = 16
ecc = 4
ecc_decimation = size // ecc
w = GF(roots_of_unity[round(math.log2(size))])
assert w**size == 1
assert (w**ecc_decimation) ** ecc == 1


pntt_shift = []
for i in range(size):
    blk = (i // ecc) * ecc
    j = i - blk
    pntt_shift.append(w ** (j * rbo(size, blk)))


def partial_ntt(buf):
    for i in range(ecc_decimation):
        buf[i * ecc : (i + 1) * ecc] = ntt(w**ecc_decimation, buf[i * ecc : (i + 1) * ecc])

    for i in range(size):
        buf[i] *= pntt_shift[i]

    for i in range(1, ecc_decimation):
        buf[:ecc] = [a + b for a, b in zip(buf[0:ecc], buf[i * ecc : (i + 1) * ecc])]

    return buf[:ecc]


def partial_ntt_block(ecc, shift):
    ecc = ntt(w**ecc_decimation, ecc)
    ecc = [e * w ** (shift * i) for i, e in enumerate(ecc)]
    return ecc


def test():
    b = [1000 + i * 10 for i in range(size)]
    print("b:", to_int_list(b))

    B = ntt(w, b)
    print("B:", to_int_list(B))
    print()

    b1 = b[:ecc] + [0] * (size - ecc)
    print("b1:", to_int_list(b1))

    b2 = [0] * ecc + b[ecc : 2 * ecc] + [0] * (size - ecc * 2)
    print("b2:", to_int_list(b2))

    b3 = [0] * ecc * 2 + b[2 * ecc : 3 * ecc] + [0] * (size - ecc * 3)
    print("b3:", to_int_list(b3))

    b4 = [0] * ecc * 3 + b[3 * ecc :]
    print("b4:", to_int_list(b4))

    print()

    B1 = ntt(w, b1)
    B11 = partial_ntt_block(b1[:ecc], shift=0)
    print("B1:", to_int_list(B1))
    print("B1:", to_int_list(B11))
    B2 = ntt(w, b2)
    B21 = partial_ntt_block(b[ecc : 2 * ecc], shift=rbo(size, ecc))
    print("B2:", to_int_list(B2))
    print("B2:", to_int_list(B21))
    B3 = ntt(w, b3)
    B31 = partial_ntt_block(b[2 * ecc : 3 * ecc], shift=rbo(size, 2 * ecc))
    print("B3:", to_int_list(B3))
    print("B3:", to_int_list(B31))
    B4 = ntt(w, b4)
    B41 = partial_ntt_block(b[3 * ecc :], shift=rbo(size, 3 * ecc))
    print("B4:", to_int_list(B4))
    print("B4:", to_int_list(B41))

    print()

    BS = [b1 + b2 + b3 + b4 for b1, b2, b3, b4 in zip(B1, B2, B3, B4)]
    print("BS:", to_int_list(BS))

    BS1 = [b1 + b2 + b3 + b4 for b1, b2, b3, b4 in zip(B11, B21, B31, B41)]
    print("BS1:", to_int_list(BS1))

    assert BS == B
    assert BS1 == B[:ecc]

    B1 = partial_ntt(b)
    print("B1:", to_int_list(B1))

    assert B1 == BS1


if __name__ == "__main__":
    test()
