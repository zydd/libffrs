import random

from ffrs.reference import *
from ffrs.reference.ntt import *
from ffrs.reference.linalg import *

GF256 = GF(2, 8)

pru = GF256.primitive_roots_of_unity()
print(pru.keys())

n = 5
arr = [GF256(x) for x in random.randbytes(n)]




w = GF256(random.choice(pru[n]))
print(f"w: {int(w)}  wi: {int(w.inv())}")


sru = [x for x in pru[n] if GF256(x).pow(2) == w]
psi = GF256(next(iter(sru)))

w_exp = [w.pow(i) for i in range(n+1)]
c = circulant(w_exp)
print_mat(c)
ci = inverse(GF256, c)
print_mat(ci)
assert identity(GF256, len(c)) == matmul(GF256, c, ci)

for i in range(1, len(c)):
    c[i] = [x * w.inv().pow(i) for x in c[i]]

print_mat(c)


for i in range(len(c)):
    c[i] = [w_exp.index(int(x)) for x in c[i]]
print("log:")
print_mat(c)

quit()

