from ffrs.reference import P
from ffrs.reference.linalg import *

def ntt(GF, w, arr):
    poly = P(GF, arr)
    return [poly.eval(w.pow(i)) for i in range(len(arr))]


def intt(GF, w, arr):
    return ntt(GF, w.inv(), arr)


def vandermonde_ntt(w, n):
    res = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            res[i][j] = w.pow(i * j)
    return res


def vandermonde_neg_ntt(p, n):
    res = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            res[i][j] = p.pow(2 * i * j + i)
    return res


def vandermonde_neg_intt(p, n):
    res = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            res[i][j] = p.pow(-2 * i * j - j)
    return res


# def intt(GF, w, arr):
#     wi = w.inv()

#     res = [None] * len(arr)

#     for j in range(len(arr)):
#         res[j] = GF(0)
#         for i in range(len(arr)):
#             res[j] += wi.pow(i * j % len(arr)) * arr[i]
#         # res[j] /= len(arr)

#     return res

