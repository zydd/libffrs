import random

from ffrs.reference import *
from ffrs.reference.linalg import *

GF256 = GF(2, 8)

pru = GF256.primitive_roots_of_unity()

n = 17
arr = [GF256(x) for x in random.randbytes(n)]

w = GF256(random.choice(pru[n]))
print(f"w: {int(w)}  wi: {int(w.inv())}")

def encode(message, ecc_len):
    message = list(message)
    codeword = [GF256(0)] * ecc_len

    for j in range(ecc_len):
        codeword[j] = GF256(0)
        for i in range(len(message)):
            ex = (j - i) % len(message)
            codeword[j] += w.pow(ex) * message[i]
            print(" ", ex, end="")
        print()

    return message + codeword


def synds(message, ecc_len):
    message_enc = encode(message[:-ecc_len], ecc_len)
    return [a - b for a, b in zip(message[-ecc_len:], message_enc[-ecc_len:])]


def synd_ref(err_pos, err_mag, msg_len, ecc_len):
    print("err_pos", err_pos)
    print("err_mag", err_mag)

    synd = [GF256(0)] * ecc_len

    for j in range(ecc_len):
        synd[j] = GF256(0)
        for i, e in zip(err_pos, err_mag):
            ex = (j - i) % msg_len
            synd[j] += w.pow(ex) * GF256(e)

    print("synd_ref", [int(x) for x in synd])
    return synd


def decode_erasure(message, synd, ecc_len, err_pos):
    Xs = [GF256.a.pow(i) for i in err_pos]
    print("Xs", [int(i.log()) for i in Xs])

    mat = [[xk.pow(i) for xk in Xs] for i in range(len(err_pos))]

    Ys = gaussian_elim(mat, synd[:len(err_pos)])
    # print([[int(col) for col in row] for row in mat])

    print("Ys", [int(i) for i in Ys])

    for i, e in zip(err_pos, Ys):
        message[i] -= e

    return message[:-ecc_len]


def find_errors(synd, n):
    ne = len(synd)//2

    mat = [synd[i:i+ne] for i in range(ne)]

    print_mat(mat)

    err_loc_coefs = gaussian_elim(mat, synd[ne:2*ne])
    lm = P(GF256, [1] + err_loc_coefs[::-1])

    roots = [i for i in range(n) if int(lm.eval(GF256.a.pow(i).inv())) == 0]
    return roots



ecc_len = 4
message_len = 15 - ecc_len

pru = GF256.primitive_roots_of_unity()
root = GF256.a
root_i = root.inv()

assert message_len + ecc_len <= 255

for _ in range(100):
    message = [GF256(x) for x in random.randbytes(message_len)]
    # print([int(x) for x in message])

    message_enc = encode(message, ecc_len)
    err_pos = set()
    while len(err_pos) < ecc_len//2:
        err_pos.add(random.randint(0, message_len + ecc_len - 1))
    err_mag = set()
    while len(err_mag) < ecc_len//2:
        err_mag.add(random.randint(1, 255))
    err_pos = list(err_pos)
    err_mag = list(err_mag)

    if max(err_pos) >= message_len:
        print("* Warning: corrupted codeword")

    synds_ref = synd_ref(err_pos, err_mag, message_len, ecc_len)

    for i, e in zip(err_pos, err_mag):
        message_enc[i] += GF256(e)

    synd = synds(message_enc, ecc_len)
    print(f"synd     {[int(x) for x in synd]}")

    assert synd == synds_ref

    err_pos_rec = find_errors(synd, len(message_enc))
    print(f"err_pos_rec {[int(x) for x in err_pos_rec]}")

    message_dec = decode_erasure(message_enc, synd, ecc_len, err_pos_rec)
    assert message == message_dec

    # print([int(x) for x in message_dec])
