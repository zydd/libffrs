# https://cryptographycaffe.sandboxaq.com/posts/ntt-02/

import collections
import itertools
import math
import random

from ffrs.reference import F
import  ffrs.reference.ntt as ref_ntt

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

n = 8
q = next(prime_root(2**16))
q = 65537
# q = 167772161  # 167772161 % 2**25 == 1

print("q:", q)

GF = lambda x: F(q, x)


# Nth root of unity
w = next(i for i in range(2, q) if i**n % q == 1)
# 2Nth root of unity
psi = next(i for i in range(2, q) if i**2 % q == w and i**n % q == -1 % q)


print("w:", w)
print("psi:", psi)
# w = psi
# assert pow(w, n, q) == 1


def naive_ntt(a, w=w, q=q):
    out = [0] * len(a)

    for i in range(len(a)):
        for j in range(len(a)):
            out[i] = (out[i] + a[j] * w ** (i * j) % q) % q

    return out

def naive_intt(a, gen=w, modulus=q):
  deg_d = len(a)
  out = [0] * deg_d

  omegas = [0] * deg_d
  omegas[0] = 1
  for i in range(1, len(omegas)):
    omegas[i] = omegas[i-1] * pow(gen, -1, modulus) % modulus

  for i in range(deg_d):
    for j in range(deg_d):
      out[i] = (out[i] + a[j] * omegas[i * j % deg_d]) % modulus

  # Scale it down before returning.
  scaler = pow(deg_d, -1, modulus)
  return [i * scaler % modulus for i in out]


def brv(x, n):
    """ Reverses a n-bit number """
    return int(''.join(reversed(bin(x)[2:].zfill(n))), 2)


def ntt_iter(a, w=w, q=q):
    deg_d = len(a)

    # Start with stride = 1.
    stride = 1

    # Shuffle the input array in bit-reversal order.
    nbits = int(math.log2(deg_d))
    res = [a[brv(i, nbits)] for i in range(deg_d)]

    # Pre-compute the generators used in different stages of the recursion.
    gens = [w ** (2**i) % q for i in range(nbits)]
    # The first layer uses the lowest (2nd) root of unity, hence the last one.
    gen_ptr = len(gens) - 1

    # Iterate until the last layer.
    while stride < deg_d:
        # For each stride, iterate over all N//(stride*2) slices.
        for start in range(0, deg_d, stride * 2):
            # For each pair of the CT butterfly operation.
            for i in range(start, start + stride):
                # Compute the omega multiplier. Here j = i - start.
                zp = gens[gen_ptr] ** (i - start) % q

                # Cooley-Tukey butterfly.
                a = res[i]
                b = res[i+stride]
                res[i] = (a + zp * b) % q
                res[i+stride] = (a - zp * b) % q

        # Grow the stride.
        stride <<= 1
        # Move to the next root of unity.
        gen_ptr -= 1

    return res


def intt_iter(a, gen=w, q=q):
    deg_d = len(a)

    # Start with stride = N/2.
    stride = deg_d // 2

    # Shuffle the input array in bit-reversal order.
    nbits = int(math.log2(deg_d))
    res = a[:]

    # Pre-compute the inverse generators used in different stages of the recursion.
    gen = pow(gen, -1, q)
    gens = [gen ** (2**i) % q for i in range(nbits)]
    # The first layer uses the highest (d-th) root of unity, hence the first one.
    gen_ptr = 0

    # Iterate until the last layer.
    while stride > 0:
            # For each stride, iterate over all N//(stride*2) slices.
            for start in range(0, deg_d, stride * 2):
                # For each pair of the CT butterfly operation.
                for i in range(start, start + stride):
                    # Compute the omega multiplier. Here j = i - start.
                    zp = gens[gen_ptr] ** (i - start) % q

                    # Gentleman-Sande butterfly.
                    a = res[i]
                    b = res[i+stride]
                    res[i] = (a + b) % q
                    res[i+stride] = ((a - b) * zp) % q

            # Grow the stride.
            stride >>= 1
            # Move to the next root of unity.
            gen_ptr += 1

    # Scale it down before returning.
    scaler = pow(deg_d, -1, q)

    # Reverse shuffle and return.
    return [(res[brv(i, nbits)] * scaler) % q for i in range(deg_d)]


a = [random.randint(0, q - 1) for _ in range(n)]
print("a:   ", a)
print()

a_ntt_naive = naive_ntt(a)
a_ntt_ref = [int(x) for x in ref_ntt.ntt(GF, GF(w), a)]
a_ntt = ntt_iter(a)
print("naive:  ", a_ntt_naive)
print("ref:    ", a_ntt_ref)
print("iter:   ", a_ntt)
print()

assert a_ntt_naive == a_ntt_ref == a_ntt

a_intt_naive = naive_intt(a_ntt_naive)
a_intt_ref = [int(x) for x in ref_ntt.intt(GF, GF(w), a_ntt)]
a_intt = intt_iter(a_ntt)
print("naive: ", a_intt_naive)
print("ref:   ", a_intt_ref)
print("iter:  ", a_intt)
print()

assert a == a_intt_naive
assert a == a_intt_ref
assert a == a_intt
