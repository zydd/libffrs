import collections
import itertools
import math
import random

import  ffrs.reference.ntt as ref_ntt
import ffrs.reference
from ffrs.reference.util import *
import ffrs.reference.ntt as ntt
import ffrs.reference.linalg as lin


def primes():
    marked = collections.defaultdict(list)
    primes = []
    for x in itertools.count(2):
        if x in marked:
            for p in marked.pop(x):
                marked[x + p].append(p)
        else:
            primes.append(x)
            marked[x * x].append(x)
            yield x


def prime_root(n):
    for p in primes():
        if p < n: continue
        # print("p:", p, end='\r')
        if p % n == 1: yield p


def naive_ntt(a, w, q):
    out = [0] * len(a)

    for i in range(len(a)):
        for j in range(len(a)):
            out[i] = (out[i] + a[j] * w ** (i * j) % q) % q

    return out


def naive_intt(a, gen, modulus):
    size = len(a)
    out = [0] * size

    omegas = [0] * size
    omegas[0] = 1
    for i in range(1, len(omegas)):
        omegas[i] = omegas[i-1] * pow(gen, -1, modulus) % modulus

    for i in range(size):
        for j in range(size):
            out[i] = (out[i] + a[j] * omegas[i * j % size]) % modulus

    # Scale it down before returning.
    scaler = pow(size, -1, modulus)
    return [i * scaler % modulus for i in out]


def ct_ntt_iter(a, root, q, start_pos=None, end_pos=None):
    """https://cryptographycaffe.sandboxaq.com/posts/ntt-02/"""
    size = len(a)

    if start_pos is None:
        start_pos = 0

    if end_pos is None:
        end_pos = size

    # Start with stride = 1.
    stride = 1

    # Shuffle the input array in bit-reversal order.
    nbits = int(math.log2(size))
    res = rbo_sorted(a)

    roots = [1]
    for i in range(1, q):
        roots.append(roots[-1] * root % q)

    # Pre-compute the generators used in different stages of the recursion.
    # gens = [root ** (2**i) % q for i in range(nbits)]
    # The first layer uses the lowest (2nd) root of unity, hence the last one.
    gen_ptr = nbits - 1
    exp_f = size // 2

    # Iterate until the last layer.
    while stride < size:
        # For each stride, iterate over all N//(stride*2) slices.
        # print("-" * 10)
        for start in filter(lambda x: x >= start_pos, range(0, end_pos, stride * 2)):
            # For each pair of the CT butterfly operation.
            # print(f"butterfly: [{start:2} {start + stride:2}]  [{start + stride:2} {start + 2 * stride:2}]")
            for i in range(start, start + stride):
                # Compute the omega multiplier. Here j = i - start.
                # zp = root ** ((i - start) << gen_ptr) % q
                zp = roots[exp_f * (i - start)]

                # Cooley-Tukey butterfly.
                a = res[i]
                b = res[i+stride]
                res[i] = (a + zp * b) % q
                res[i+stride] = (a - zp * b) % q

        # Grow the stride.
        stride <<= 1
        # Move to the next root of unity.
        gen_ptr -= 1
        exp_f //= 2

    return res


def gs_ntt_iter(a, root, q, end_pos=None):
    """https://cryptographycaffe.sandboxaq.com/posts/ntt-02/"""
    size = len(a)

    if end_pos is None:
        end_pos = size

    res = a[:]
    exp_f = 0

    roots = [1]
    for i in range(1, q):
        roots.append(roots[-1] * root % q)

    stride = size // 2
    while stride > 0:
        # For each stride, iterate over all N//(stride*2) slices.
        print("-" * 10)
        for start in range(0, size, stride * 2):
            # For each pair of the CT butterfly operation.
            if start >= end_pos:
                continue
            print(f"butterfly: [{start:2} {start + stride:2}]")
            for i in range(start, start + stride):
                # Compute the omega multiplier. Here j = i - start.
                # zp = root ** ((i - start) << exp_f) % q
                zp = roots[(i - start) << exp_f]
                
                # Gentleman-Sande butterfly.
                a = res[i]
                b = res[i+stride]
                res[i] = (a + b) % q
                res[i+stride] = ((a - b) * zp) % q

        # Grow the stride.
        stride >>= 1
        # Move to the next root of unity.
        exp_f += 1

    return rbo_sorted(res)


n = 16
# q = next(prime_root(2**16))
q = 65537

print("q:", q)

GF = ffrs.reference.GF(q)


# Nth root of unity
w = next(i for i in range(2, q) if i**n % q == 1)
# 2Nth root of unity
psi = next(i for i in range(2, q) if i**2 % q == w and i**n % q == -1 % q)


print("w:", w)
print("psi:", psi)
w = int(GF(w).inv())

# test data
a = [random.randint(0, q - 1) for _ in range(n)]
print("a:   ", a)
print()

a_ntt_naive = naive_ntt(a, w, q)
a_ntt_ref = to_int_list(ref_ntt.ntt(GF, GF(w), a))
a_ntt = ct_ntt_iter(a, w, q)
a_ntt2 = gs_ntt_iter(a, w, q)
print("naive:  ", a_ntt_naive)
print("ref:    ", a_ntt_ref)
print("ct-iter:", a_ntt)
print("gs-iter:", a_ntt2)

v = ntt.vandermonde_ntt(GF(w), len(a))
v_iref = ntt.vandermonde_intt(GF(w).inv(), GF(len(a)))
v_i = lin.inverse(GF, v)
assert v_i == v_iref

a_ntt_v = lin.matmul(GF, v, GF(a))
print("ntt v:  ", to_int_list(a_ntt_v))

print()


a_intt_naive = naive_intt(a_ntt, w, q)
a_intt_ref = to_int_list(ref_ntt.intt(GF, GF(w), a_ntt))
a_intt_ct = to_int_list([GF(x) // GF(len(a)) for x in ct_ntt_iter(a_ntt, int(GF(w).inv()), q)])
a_intt_gs = to_int_list([GF(x) // GF(len(a)) for x in gs_ntt_iter(a_ntt, int(GF(w).inv()), q)])
a_intt_v = to_int_list(lin.matmul(GF, v_i, GF(a_ntt)))
print("naive:  ", a_intt_naive)
print("ref:    ", a_intt_ref)
print("ct-iter:", a_intt_ct)
print("gs-iter:", a_intt_gs)
print("intt v: ", a_intt_v)
print()


assert a_ntt_naive == a_ntt_ref == a_ntt == a_ntt_v
assert a == a_intt_naive == a_intt_ref ==  a_intt_ct == a_intt_v == a_intt_gs



end = n//2
gs_partial = gs_ntt_iter(a, w, q, end_pos=end)
ct_partial = ct_ntt_iter(a, w, q)
print("a:   ", gs_partial)
print("b:   ", ct_partial)
assert gs_partial[:end] == ct_partial[:end]
