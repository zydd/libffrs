from ffrs.reference import *
from ffrs.reference.ntt import *
from ffrs.reference.linalg import *

GF2 = GF(2)

m = P(GF2, 0xabcdef12)

poly = P(GF2, 0x11d)
print("ref div", hex(GF.poly_to_int(2, m // poly)))
print("ref rem", hex(GF.poly_to_int(2, m % poly)))

assert (m // poly) * poly + (m % poly) == m


t = P(GF2, 1 << 32)
t_p = t // poly

div = m * t_p // t
print("m", hex(GF.poly_to_int(2, m)))
print("(t/p)", hex(GF.poly_to_int(2, t_p)))
print("m * (t/p)", hex(GF.poly_to_int(2, m * t_p)))

print("div", hex(GF.poly_to_int(2, div)))
assert div == m // poly

rem = m - div * poly
print("rem", hex(GF.poly_to_int(2, rem)))
assert rem == m % poly
0xcddd2dae724d58eb50f9d5f0