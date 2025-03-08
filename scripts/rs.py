import random

from ffrs.reference import *
from ffrs.reference.linalg import *

GF256 = GF(2, 8)


def encode(message, ecc_len):
    message = list(message)
    message_poly = P(GF256, message)
    for i in range(ecc_len):
        message.append(message_poly.eval(GF256.a.pow(i)))
    return message


def encodev2(message, ecc_len):
    message = list(message)
    message_poly = P(GF256, message)
    evaluated_message = []

    for i in range(ecc_len):
        evaluated_message.append(message_poly.eval(root.pow(i)))

    for i in range(ecc_len):
        assert root_i.pow(i * len(message)) == root_i.pow(i).pow(len(message))
        evaluated_message[i] *= root_i.pow(i * len(message))

    codeword = matmul(GF256, v_i, evaluated_message)

    return message + codeword


def syndsv2(message_enc, ecc_len):
    message = message_enc[:-ecc_len]
    codeword = message_enc[-ecc_len:]
    message_poly = P(GF256, message)
    evaluated_message = []

    for i in range(ecc_len):
        evaluated_message.append(message_poly.eval(root.pow(i)))

    codeword_rec = matmul(GF256, v, codeword)

    for i in range(ecc_len):
        codeword_rec[i] *= root.pow(i * len(message))

    synd = vec_add(evaluated_message, codeword_rec)

    return synd


def synds(message, ecc_len):
    message_poly = P(GF256, message[:-ecc_len])
    return [message_poly.eval(GF256.a.pow(i)) - message[-ecc_len + i]
                for i in range(ecc_len)]


def synd_ref(err_pos, err_mag, ecc_len):
    print("err_pos", err_pos)
    print("err_mag", err_mag)

    synd = []

    for i in range(ecc_len):
        synd.append(GF256(0))
        for pos, mag in zip(err_pos, err_mag):
            synd[-1] += GF256.a.pow(i * pos) * GF256(mag)

    print("synd_ref", [int(x) for x in synd])


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

    err_loc_coefs = gaussian_elim(mat, synd[ne:2*ne])
    lm = P(GF256, [1] + err_loc_coefs[::-1])

    roots = [i for i in range(n) if int(lm.eval(GF256.a.pow(i).inv())) == 0]
    return roots



ecc_len = 4
message_len = 15 - ecc_len

pru = GF.primitive_roots_of_unity(GF256)
root = GF256.a
root_i = root.inv()
v = vandermonde_ntt(root, ecc_len)
# v_i = vandermonde_ntt(root_i, ecc_len)
v_i = inverse(GF256, v)

assert matmul(GF256, v, v_i) == identity(GF256, ecc_len)

assert message_len + ecc_len <= 255

for _ in range(100):
    message = [GF256(x) for x in random.randbytes(message_len)]
    # print([int(x) for x in message])

    message_enc = encodev2(message, ecc_len)
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

    synd_ref(err_pos, err_mag, ecc_len)

    for i, e in zip(err_pos, err_mag):
        message_enc[i] += GF256(e)

    synd = syndsv2(message_enc, ecc_len)
    print(f"synd     {[int(x) for x in synd]}")

    err_pos_rec = find_errors(synd, len(message_enc))
    print(f"err_pos_rec {[int(x) for x in err_pos_rec]}")

    message_dec = decode_erasure(message_enc, synd, ecc_len, err_pos_rec)
    assert message == message_dec

    # print([int(x) for x in message_dec])
