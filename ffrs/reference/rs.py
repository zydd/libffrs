import ffrs.reference
import ffrs.reference.ntt

from ffrs.reference import P
from ffrs.reference.util import rbo, rbo_sorted, to_int_list

_intt = lambda GF, w, x: rbo_sorted(to_int_list(ffrs.reference.ntt.intt(GF, w, x)))


def locator(GF, w, err_pos):
    # print("errors:", err_pos)
    # print("root:", w)
    poly = P(GF, [1])
    for pos in err_pos:
        x = GF(w**pos)
        # print("pos", pos)
        # print("loc", poly)
        poly = poly * P(GF, [1, -x])
    return poly


def poly_roots(GF, w, poly):
    w = GF(w)
    roots = []
    for i in range(GF.field_elements):
        if poly.eval(w**-i) == GF(0):
            roots.append(i)
    return roots


def sugiyama(GF, synds):
    R2 = P(GF, [0, 1]) ** len(synds)
    R1 = P(GF, (synds))
    A2 = P(GF, [0])
    A1 = P(GF, [1])
    while R1.deg() >= len(synds) // 2:
        Q = R2 // R1

        t = A2 - Q * A1
        A2 = A1
        A1 = t

        t = R2 - Q * R1
        R2 = R1
        R1 = t

        print("R1", list(map(int, R1.x)))
        print("A1", list(map(int, A1.x)))

    locator = P(GF, [a // GF(A1.x[0]) for a in A1.x])
    evaluator = P(GF, [a // GF(A1.x[0]) for a in R1.x])
    print("locator", locator)
    print("evaluator", evaluator)
    return locator, evaluator


def forney(size, root, locator, evaluator):
    err_pos = [i for i in range(size) if locator.eval(root**-i) == 0]

    err_val = []
    for pos in err_pos:
        err_val.append(-(root**pos) * evaluator.eval(root**-pos) // locator.deriv().eval(root**-pos))

    err_pos = [rbo(size, i) for i in err_pos]
    print("err_pos:", err_pos)
    print("err_val:", list(map(int, err_val)))
    return err_pos, err_val


def mix_ecc(rs, GF, ecc):
    i = rbo(rs.block_len, rs.block_len - rs.ecc_len)
    w_i = rs.gf.pow(rs.gf.inv(rs.root), i)

    ecc_mix = [rs.gf.mul(s, rs.gf.sub(0, rs.gf.pow(w_i, j))) for j, s in enumerate(ecc)]

    # ecc_root = GF(rs.root) ** (rs.block_size // rs.ecc_size)
    ecc_root = GF(rs.gf.exp(rs.gf.div(rs.gf.log(1), rs.ecc_len)))
    ecc_mix = _intt(GF, ecc_root, ecc_mix)

    return ecc_mix
