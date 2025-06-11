from ffrs.reference import P
from ffrs.reference.linalg import *

def ntt(GF, w, arr):
    poly = P(GF, arr)
    return [poly.eval(w ** i) for i in range(len(arr))]


def intt(GF, w, arr):
    return [x // GF(len(arr)) for x in ntt(GF, w.inv(), arr)]


def vandermonde_ntt(w, n):
    res = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            res[i][j] = w ** (i * j)
    return res


def vandermonde_intt(w, n):
    n_int = int(n)
    res = [[None for _ in range(n_int)] for _ in range(n_int)]
    for i in range(n_int):
        for j in range(n_int):
            res[i][j] = w ** (i * j) // n
    return res


def vandermonde_neg_ntt(p, n):
    res = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            res[i][j] = p ** (2 * i * j + i)
    return res


def vandermonde_neg_intt(p, n):
    res = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            res[i][j] = p ** (-2 * i * j - j)
    return res


# def intt(GF, w, arr):
#     wi = w.inv()

#     res = [None] * len(arr)

#     for j in range(len(arr)):
#         res[j] = GF(0)
#         for i in range(len(arr)):
#             res[j] += wi ** (i * j % len(arr)) * arr[i]
#         # res[j] /= len(arr)

#     return res

