import math
import random
import ffrs.reference as ref
import ffrs.reference.ntt as ref_ntt
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted


GF = ref.GF(0x10001, 1, 3)
ntt = lambda w, x: rbo_sorted(ref_ntt.ntt(GF, w, x))
intt = lambda w, x: rbo_sorted(ref_ntt.intt(GF, w, x))

ntt_rbo = lambda w, x: ref_ntt.ntt(GF, w, rbo_sorted(x))
ntt_norm = lambda w, x: ref_ntt.ntt(GF, w, x)


roots_of_unity = [1, 1, 256, 4096, 64, 8, 8224, 13987, 282, 15028, 19139, 3668, 11088, 6561, 81, 9]
# roots_of_unity ** 2 = [1, 1, 1, 256, 4096, 64, 8, 8224, 13987, 282, 15028, 19139, 3668, 11088, 6561, 81]


size = 16
ecc = 4
ecc_decimation = size // ecc
w = GF(roots_of_unity[round(math.log2(size))])
assert w**size == 1

def _mix_ecc(ecc_mix):
    # i = 1
    # i = size - ecc
    # i = rbo(size, i)

    ecc_orig = ecc_mix[:]
    print("ecc:", ecc_orig)

    # ecc_mix = [GF(s) * (-w**(-ecc_decimation * i * j)) for j, s in enumerate(ecc_mix)]

    ecc_mix = ecc_mix # * ecc_decimation
    ecc_mix = intt(w ** ecc_decimation, ecc_mix)

    # ecc_mix = [-s for s in ecc_mix]
    # print("ecc_mix:", ecc_mix)
    ecc_mix = ecc_mix[:ecc]

    assert ecc_orig == _unmix_ecc(ecc_mix), "ECC mixing failed"

    return ecc_mix

def _unmix_ecc(ecc_mix):
    ecc_unmix = ecc_mix
    # ecc_unmix = [0] * (size - ecc) + ecc_unmix

    # i = size - ecc
    # i = rbo(size, i)
    # ecc_unmix = [GF(s) * (w**(ecc_decimation * i * j)) for j, s in enumerate(ecc_unmix)]

    ecc_unmix = ntt_rbo(w ** ecc_decimation, ecc_unmix)
    ecc_unmix = ecc_unmix[:ecc]
    # ecc_unmix = rbo_sorted(ecc_unmix)

    # i = size - ecc
    # i = rbo(ecc, i)
    print("ecc_unmix:", ecc_unmix)
    return ecc_unmix


def encode(buf):
    buf = GF(buf)
    assert w**size == 1
    ecc_mix = ntt(w, buf)[:ecc]

    # ecc_mix = _mix_ecc(ecc_mix)

    buf[-ecc:] = ecc_mix

    return buf


def find_errors(msg1):
    msg2 = msg1[:-ecc] + [0] * ecc
    ecc1 = msg1[-ecc:]
    ecc2 = encode(msg2)[-ecc:]

    synds = [a - b for a, b in zip(GF(ecc1), GF(ecc2))]
    # synds = _unmix_ecc(synds)
    # synds = ntt(w, msg1)[:ecc]
    synds = rbo_sorted(synds)
    print("synds:", synds)

    if all(s == 0 for s in synds):
        return {}

    for err_count in range(ecc // 2, 0, -1):
        try:
            mat = [synds[i:i+err_count] for i in range(err_count)]

            err_loc_coefs = ref.linalg.gaussian_elim(mat, [-s for s in synds[err_count:2*err_count]])
            lm = ref.P(GF, [1] + err_loc_coefs[::-1])

            err_pos = [i for i in range(size * 2) if int(lm.eval(w ** -i)) == 0]
            if err_pos:
                break
        except:
            continue
    else:
        raise RuntimeError("Decoding failed")


    err_ws = [[w ** (err_pos[i] * j) for i in range(err_count)] for j in range(err_count)]
    # err_pos = [rbo(size, i) for i in err_pos]
    err_pos = [i // ecc_decimation for i in err_pos]

    err_mag = ref.linalg.gaussian_elim(err_ws, GF(synds[:err_count]))
    errors = dict(zip(err_pos, map(int, err_mag)))
    return errors


def test():
    buf = [random.randrange(0, 2**16) for _ in range(size - ecc)] + [0] * ecc
    print("buf:", buf)
    print()

    res = encode(buf)

    res[2] -= GF(444)
    res[3] -= GF(666)
    # res[-ecc + 0] += GF(123)
    # res[-ecc + 3] += GF(888)

    errs = find_errors(res)
    print()
    print("errors:", errs)


if __name__ == "__main__":
    test()
