#
#  test_lib_gfi16_ntt.py
#
#  Copyright 2026 Gabriel Machado
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import math
import pytest

import ffrs

import ffrs.reference as ref
from ffrs.reference.util import randbytes, to_int_list

GF = ref.GF(65537, 1, 3)
gf = ffrs.GFi16(3)


P = lambda c: ref.P(GF, c)


@pytest.mark.parametrize(
    "ntt",
    [
        ffrs.NTT(gf, 4, 2),
        ffrs.NTT(gf, 8, 4),
        ffrs.NTT(gf, 16, 8),
        # ffrs.NTT(gf, 32, 16),
    ],
)
class TestNTT:
    def test_inv(self, ntt: ffrs.NTT):
        data = randbytes(ntt.ntt_size)
        data_ntt = ntt.ntt(data)
        data_intt = ntt.intt(data_ntt)
        assert data != data_ntt
        assert data == data_intt

    def test_poly_mul(self, ntt: ffrs.NTT):
        poly_a = randbytes(ntt.ecc_size)
        poly_b = randbytes(ntt.ecc_size)
        poly_a[ntt.ecc_size // 2 :] = bytearray(ntt.ecc_size // 2)
        poly_b[ntt.ecc_size // 2 :] = bytearray(ntt.ecc_size // 2)

        res = ntt.poly_mul(poly_a, poly_b)

        poly_a_ref = P(to_int_list(poly_a))
        poly_b_ref = P(to_int_list(poly_b))
        res_ref = poly_a_ref * poly_b_ref

        assert P(to_int_list(res)) == res_ref

    @pytest.mark.parametrize("n", [4, 8, 16])
    def test_poly_mul_n(self, ntt: ffrs.NTT, n):
        poly_a = randbytes(ntt.ecc_size * n)
        poly_b = randbytes(ntt.ecc_size * n)

        for i in range(n):
            start = i * ntt.ecc_size
            poly_a[start + ntt.ecc_size // 2 : start + ntt.ecc_size] = bytearray(ntt.ecc_size // 2)
            poly_b[start + ntt.ecc_size // 2 : start + ntt.ecc_size] = bytearray(ntt.ecc_size // 2)

        poly_mul = getattr(ntt, f"poly_mul{n}")
        res = poly_mul(poly_a, poly_b)

        for i in range(n):
            start = i * ntt.ecc_size
            res_ref = ntt.poly_mul(poly_a[start : start + ntt.ecc_size], poly_b[start : start + ntt.ecc_size])
            assert res[start : start + ntt.ecc_size] == res_ref

    def test_poly_div_zero(self, ntt: ffrs.NTT):
        poly_a = randbytes(ntt.ecc_size)

        # A // 0
        with pytest.raises(RuntimeError):
            res = ntt.poly_div(poly_a, bytearray(ntt.ecc_size))

        # 0 // 0
        with pytest.raises(RuntimeError):
            res = ntt.poly_div(bytearray(ntt.ecc_size), bytearray(ntt.ecc_size))

    @pytest.mark.skip
    def test_poly_div_same_size(self, ntt: ffrs.NTT, subtests):
        for size in range(2, ntt.ecc_size + 1, 2):
            with subtests.test(size=size):
                poly_a = bytearray(ntt.ecc_size)
                poly_b = bytearray(ntt.ecc_size)

                poly_a[:size] = randbytes(size)
                poly_b[:size] = randbytes(size)

                res = ntt.poly_div(poly_a, poly_b)

                poly_a_ref = P(to_int_list(poly_a))
                poly_b_ref = P(to_int_list(poly_b))
                res_ref = poly_a_ref // poly_b_ref

                assert P(to_int_list(res)) == res_ref

    def test_poly_div_higher_deg_denominator(self, ntt: ffrs.NTT, subtests):
        for size in range(2, ntt.ecc_size - 1, 2):
            with subtests.test(size=size):
                poly_a = bytearray(ntt.ecc_size)
                poly_b = bytearray(ntt.ecc_size)

                poly_a[:size] = randbytes(size)
                poly_b[: size + 2] = randbytes(size + 2)

                res = ntt.poly_div(poly_a, poly_b)

                assert P(to_int_list(res)) == P([0])

    @pytest.mark.skip
    def test_poly_div_lower_deg_denominator(self, ntt: ffrs.NTT, subtests):
        for num_size in range(2, ntt.ecc_size + 1, 2):
            for den_size in range(2, num_size, 2):
                with subtests.test(num_size=num_size, den_size=den_size):
                    poly_a = bytearray(ntt.ecc_size)
                    poly_b = bytearray(ntt.ecc_size)

                    poly_a[:num_size] = randbytes(num_size)
                    poly_b[:den_size] = randbytes(den_size)

                    res = ntt.poly_div(poly_a, poly_b)

                    poly_a_ref = P(to_int_list(poly_a))
                    poly_b_ref = P(to_int_list(poly_b))
                    res_ref = poly_a_ref // poly_b_ref
                    res_ref.x[0] = GF(0)

                    assert P(to_int_list(res)) == res_ref

    @pytest.mark.parametrize("n", [4, 8, 16])
    def test_poly_div_n(self, ntt: ffrs.NTT, n):
        poly_a = randbytes(ntt.ecc_size * n)
        poly_b = randbytes(ntt.ecc_size * n)

        for i in range(n):
            # mirror single-block zeros used in test_poly_div
            start = i * ntt.ecc_size
            # poly_a[start + ntt.ecc_size - 2 : start + ntt.ecc_size] = bytearray(2)
            poly_b[start + ntt.ecc_size - 2 : start + ntt.ecc_size] = bytearray(2)

        poly_div = getattr(ntt, f"poly_div{n}")
        res = poly_div(poly_a, poly_b)

        for i in range(n):
            start = i * ntt.ecc_size
            res_ref = ntt.poly_div(poly_a[start : start + ntt.ecc_size], poly_b[start : start + ntt.ecc_size])
            assert res[start : start + ntt.ecc_size] == res_ref

    def test_poly_inv_mod(self, ntt: ffrs.NTT, subtests):
        for mod in range(1, ntt.ecc_len + 1):
            for poly_size in range(2, ntt.ecc_size + 1, 2):
                mod_exp = math.ceil(math.log2(mod))
                # Detect aliasing in second-to-last iteration
                aliasing = 2 ** (mod_exp - 1) + poly_size // 2 - 2 > ntt.ecc_len

                with subtests.test(mod=mod, aliasing=aliasing, poly_len=poly_size // 2):
                    poly_a = bytearray(ntt.ecc_size)
                    poly_a[:poly_size] = randbytes(poly_size)

                    if not aliasing:
                        poly_a_inv = ntt.poly_inv(poly_a, mod)
                        # A * inv(A) % mod
                        mul = P(to_int_list(poly_a)) * P(to_int_list(poly_a_inv)) % P([0] * mod + [1])
                        assert mul == P([1])
                    else:
                        with pytest.raises(RuntimeError):
                            poly_a_inv = ntt.poly_inv(poly_a, mod)

    @pytest.mark.parametrize("n", [4, 8, 16])
    def test_poly_inv_n(self, ntt: ffrs.NTT, n, subtests):
        mod = 1
        while mod < ntt.ecc_len // 2:
            with subtests.test(n=n, mod=mod):
                poly_a = randbytes(ntt.ecc_size * n)

                for i in range(n):
                    start = i * ntt.ecc_size
                    poly_a[start + ntt.ecc_size // 2 : start + ntt.ecc_size] = bytearray(ntt.ecc_size // 2)

                poly_inv = getattr(ntt, f"poly_inv{n}")
                res = poly_inv(poly_a, mod)

                for i in range(n):
                    start = i * ntt.ecc_size
                    res_ref = ntt.poly_inv(poly_a[start : start + ntt.ecc_size], mod)
                    assert res[start : start + ntt.ecc_size] == res_ref
            mod *= 2

    # TODO: mul, div, inv tests with polynomials of different degrees + simd
