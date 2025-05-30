import random

import ffrs.reference

def sample_field(GF, start=0, samples=1000):
    if GF.field_elements < samples:
        return range(GF.field_elements)

    n = 8

    # Yield first n elements
    for i in range(start, min(n, GF.field_elements)):
        yield i

    # Yield last n elements
    for i in range(max(n, GF.field_elements - n), GF.field_elements):
        yield i

    for _ in range(samples - 2 * n):
        yield random.randrange(start, GF.field_elements)


def to_bytearray(v, sizeof, byteorder="little"):
    if type(v) in [bytearray, bytes, str]:
        return bytearray(v, "utf-8")

    buf = bytearray(len(v) * sizeof)
    for i, v in enumerate(v):
        buf[i * sizeof:(i + 1) * sizeof] = int.to_bytes(v, sizeof, byteorder)
    return buf


def to_int_list(v, sizeof=None, byteorder="little"):
    if type(v) is list:
        if len(v) == 0: return []

        if type(v[0]) in [ffrs.reference.F, ffrs.reference.Fp]:
            return list(map(int, v))
        else:
            raise TypeError
    else:
        assert sizeof, "Must specify sizeof for byte objects"
        return [int.from_bytes(v[i:i+sizeof], byteorder) for i in range(0, len(v), sizeof)]



def randbytes(n, start=0, stop=256):
    return bytearray(random.randrange(start=start, stop=stop) for _ in range(n))

