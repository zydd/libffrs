import random

import ffrs
from ffrs.reference import GF
from ffrs.reference.ntt import ntt

NTT = ffrs.NTT256x8()


GF256 = GF(2, 8, NTT.gf.primitive, NTT.gf.poly1 | 0x100)


random.seed(42)

def randbytes(n, start=0, stop=256):
    return bytearray(random.randrange(start=start, stop=stop) for _ in range(n))


def test_ntt8():
    a = randbytes(8)
    res = NTT.ntt8(a)
    ref = ntt(GF256, GF256(NTT.gf.primitive), a)
    print(list(a))
    print("ref", [int(x) for x in ref])
    print("res", list(res))

    ref2 = ffrs.GF256().poly_eval8(a[::-1], bytearray(int(GF256(NTT.gf.primitive).pow(i)) for i in range(len(a))))
    print("ref2", list(ref2))
    assert [int(x) for x in ref] == list(res) == list(ref2)


if __name__ == "__main__":
    test_ntt8()
