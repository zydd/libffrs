import math
import random
import ffrs.reference as ref
import ffrs.reference.ntt as ref_ntt
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted

GF = ref.GF(0x10001, 1, 3)
P = lambda x: ref.P(GF, x)
x = P([0, 1])

nttr = lambda w, x: rbo_sorted(ref_ntt.ntt(GF, w, x))
intt = lambda w, x: rbo_sorted(ref_ntt.intt(GF, w, x))


roots_of_unity = [1, 1, 256, 4096, 64, 8, 8224, 13987, 282, 15028, 19139, 3668, 11088, 6561, 81, 9]
# roots_of_unity ** 2 = [1, 1, 1, 256, 4096, 64, 8, 8224, 13987, 282, 15028, 19139, 3668, 11088, 6561, 81]


size = 16
ecc = 4
ecc_decimation = size // ecc
w = GF(roots_of_unity[round(math.log2(size))])
assert w**size == 1
assert (w**ecc_decimation) ** ecc == 1


pintt_shift = []
for i in range(size):
    blk = (i // ecc) * ecc
    j = i - blk
    pintt_shift.append(w ** (-j * rbo(size, blk)) // GF(size))


def partial_intt_block(block, shift):
    block = [c * e for c, e in zip(pintt_shift[shift * ecc :], block)]
    block = nttr(w**-ecc_decimation, block)
    return list((block))


def test():
    b = (P([1]) - x * w**8) * (P([1]) - x * w**10)
    print(b)
    print(f"b(w**-8) =", b(w**-8))
    print(f"b(w**-10) =", b(w**-10))

    b = b.x
    b += [GF(0)] * (size - len(b))
    print("b:", to_int_list(b))

    B = intt(w, b)
    print("B:", to_int_list(B))
    print()

    print("zeros:", [rbo(size, i) for i, x in enumerate(B) if x == 0])
    print()

    b1 = b[:ecc]
    BN = []
    for i in range(0, ecc_decimation):
        B1 = partial_intt_block(b1, shift=i)
        BN.extend(B1)
        print(f"B{i}:", to_int_list(B1))

    print()
    print(to_int_list(BN))

    if BN == B:
        print("OK")
    else:
        print("\n** FAIL **")


if __name__ == "__main__":
    test()
