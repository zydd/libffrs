import random

from ffrs.reference import *
from ffrs.reference.ntt import *
from ffrs.reference.linalg import *

GF256 = GF(2, 8)

pru = GF.primitive_roots_of_unity(GF256)
print(pru.keys())


n = 5
arr = [GF256(x) for x in random.randbytes(n)]


w = GF256(random.choice(pru[n]))
w = GF256.a
print(f"w: {int(w)}  wi: {int(w.inv())}")

v = vandermonde_ntt(w, n)
v_i = inverse(GF256, v)
v_w_i = vandermonde_ntt(w.inv(), n)

arr_n = ntt(GF256, w, arr)
arr_v = matmul(GF256, v, arr)
print(f"arr  {[int(x) for x in arr]}")
print(f"ntt  {[int(x) for x in arr_n]}")
print(f"vntt {[int(x) for x in arr_n]}")

assert arr_n == arr_v

arr_r = intt(GF256, w, arr_n)
arr_vr = matmul(GF256, v_i, arr_v)
arr_vwr = matmul(GF256, v_w_i, arr_v)
print(f"intt {[int(x) for x in arr_r]}")
print(f"vintt{[int(x) for x in arr_vr]}")
print(f"wintt{[int(x) for x in arr_vwr]}")

assert arr_vwr == arr
assert arr_vr == arr
assert arr_r == arr

