from ffrs.reference import P

def ntt(GF, w, arr):
    poly = P(GF, arr)
    return [poly.eval(w.pow(i)) for i in range(len(arr))]


def intt(GF, w, arr):
    return ntt(GF, w.inv(), arr)


# def intt(GF, w, arr):
#     wi = w.inv()

#     res = [None] * len(arr)

#     for j in range(len(arr)):
#         res[j] = GF(0)
#         for i in range(len(arr)):
#             res[j] += wi.pow(i * j % len(arr)) * arr[i]
#         # res[j] /= len(arr)

#     return res

