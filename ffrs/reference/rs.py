from ffrs.reference import P
from ffrs.reference.util import rbo


def locator(GF, w, err_pos):
    poly = P(GF, [1])
    for pos in err_pos:
        x = -w ** pos
        poly = poly * P(GF, [1, x])
    return poly


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
        err_val.append(- root**pos * evaluator.eval(root**-pos) // locator.deriv().eval(root**-pos))

    err_pos = [rbo(size, i) for i in err_pos]
    print("err_pos:", err_pos)
    print("err_val:", list(map(int, err_val)))
    return err_pos, err_val

