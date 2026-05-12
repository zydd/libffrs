/**************************************************************************
 * rsi16md_impl.hpp
 *
 * Copyright 2025 Gabriel Machado
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 **************************************************************************/

# pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

#include "pygfi16.hpp"


template<typename GFT>
class RSi16vImpl {
public:
    uint32_t root;

    inline RSi16vImpl(GFi16 const& gf, size_t block_len, size_t ecc_len):
        gf(gf),
        block_len(block_len),
        ecc_len(ecc_len)
    {
        _rbo_shift = 16 - __builtin_ctzl(block_len);
        block_len_iv = GFT{0} + gf.inv(block_len);
        ecc_len_iv = GFT{0} + gf.inv(ecc_len);

        root = gf.exp(gf.div(gf.log(1), block_len));
        if (root >= 0x8000)
            root = gf.neg(root);

        if (gf.pow(root, block_len) != 1)
            throw std::runtime_error("Root of unity not found for block size");

        uint32_t root_i = gf.inv(root);
        _roots_v_block.resize(block_len);
        _roots_iv_block.resize(block_len);
        for (size_t i = 0; i < block_len; ++i) {
            _roots_v_block[i] = GFT{0} + gf.pow(root, i);
            _roots_iv_block[i] = GFT{0} + gf.pow(root_i, i);
        }

        _ecc_mix_v.resize(ecc_len);
        _ecc_mix_iv.resize(ecc_len);
        auto ecc_mix_w = *reinterpret_cast<uint32_t *>(&_roots_v_block[rbo(block_len - ecc_len)]);
        auto ecc_mix_w_i = gf.inv(ecc_mix_w);
        for (size_t i = 0; i < ecc_len; ++i) {
            _ecc_mix_v[i] = GFT{0} + gf.neg(gf.div(gf.pow(ecc_mix_w_i, i), ecc_len));
            _ecc_mix_iv[i] = GFT{0} + gf.neg(gf.pow(ecc_mix_w, i));
        }

        // uint32_t ecc_root = gf.pow(root, block_len / ecc_len);
        uint32_t ecc_root = gf.exp(gf.div(gf.log(1), ecc_len));
        uint32_t ecc_root_i = gf.inv(ecc_root);
        _roots_v_ecc.resize(ecc_len);
        _roots_iv_ecc.resize(ecc_len);
        for (size_t i = 0; i < ecc_len; ++i) {
            _roots_v_ecc[i] = GFT{0} + gf.pow(ecc_root, i);
            _roots_iv_ecc[i] = GFT{0} + gf.pow(ecc_root_i, i);
        }
    }

    inline void encode(GFT block[]) const {
        ntt(block);

        for (size_t j = 0; j < ecc_len; ++j)
            block[j] = gf.mul(block[j], _ecc_mix_v[j]);

        gs_butterfly(&_roots_iv_ecc[0], &block[0], ecc_len, ecc_len);
    }

    inline void ecc_mix(GFT ecc[]) const {
        for (size_t j = 0; j < ecc_len; ++j)
            ecc[j] = gf.mul(ecc[j], _ecc_mix_v[j]);

        gs_butterfly(&_roots_iv_ecc[0], &ecc[0], ecc_len, ecc_len);
    }

    inline void ecc_unmix(GFT ecc[]) const {
        ct_butterfly(&_roots_v_ecc[0], &ecc[0], ecc_len);

        for (size_t j = 0; j < ecc_len; ++j)
            ecc[j] = gf.mul(ecc[j], _ecc_mix_iv[j]);
    }

    inline void ntt(GFT block[]) const {
        ct_butterfly(&_roots_v_block[0], &block[0], block_len);
    }

    inline void intt(GFT block[]) const {
        gs_butterfly(&_roots_iv_block[0], &block[0], block_len, block_len);
        for (size_t j = 0; j < block_len; ++j)
            block[j] = gf.mul(block[j], block_len_iv);
    }

    inline void repair(GFT block[], GFT temp6[]) const {
        ntt(block);

        auto locator_poly = &temp6[0];
        auto locator_poly_len = sugiyama(block, locator_poly, &temp6[block_len]);
        repair_loc(block, &locator_poly[0], locator_poly_len - 1, &temp6[block_len]);

        intt(block);
    }

    inline void repair(GFT block[], const size_t error_pos_rbo[], size_t error_count, GFT temp2[]) const {
        ntt(block);

        auto locator_poly = &temp2[block_len];
        error_locator(error_pos_rbo, error_count, locator_poly);
        repair_loc(block, locator_poly, error_count, temp2);

        intt(block);
    }

    inline void repair_loc(GFT block_ntt[], const GFT locator_poly[], size_t error_count, GFT temp[]) const {
        // copy synds
        std::copy_n(&block_ntt[0], ecc_len, &temp[0]);

        for (size_t j = ecc_len; j < block_len; ++j) {
            GFT sum = GFT{0};
            for (size_t i = 0; i < error_count; ++i)
                sum = gf.sub(sum, gf.mul(locator_poly[i], temp[j - error_count + i]));

            temp[j] = sum;
            block_ntt[j] = gf.sub(block_ntt[j], sum);
        }

        std::fill_n(&block_ntt[0], ecc_len, GFT{0});

        // TODO limit to error positions only, test if performance improvement
        // for (size_t j = 0; j < error_count; ++j) {
        //     size_t i = error_pos_rbo[j];
        //     block[i] = gf.div(block[i], block_len);
        // }
    }

    uint16_t rbo(uint16_t b) const {
        return ffrs::detail::rbo16(b) >> _rbo_shift;
    }

protected:
    GFi16 const& gf;
    uint16_t _rbo_shift;
    size_t block_len;
    size_t ecc_len;
    GFT block_len_iv;
    GFT ecc_len_iv;
    std::vector<GFT> _roots_v_block;
    std::vector<GFT> _roots_iv_block;
    std::vector<GFT> _roots_v_ecc;
    std::vector<GFT> _roots_iv_ecc;
    std::vector<GFT> _ecc_mix_v;
    std::vector<GFT> _ecc_mix_iv;

    inline void ct_butterfly(const GFT roots[], GFT block[], size_t block_len) const {
        for (size_t stride = 1, exp_f = block_len >> 1; stride < block_len; stride *= 2, exp_f >>= 1) {
            for (size_t start = 0; start < block_len /*input_size*/; start += stride * 2) {
                {
                    // Cooley-Tukey butterfly
                    GFT a = block[start];
                    GFT b = block[start + stride];
                    block[start] = gf.add(a, b);
                    block[start + stride] = gf.sub(a, b);
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // j = i - start
                    GFT w = roots[exp_f * (i - start)];

                    // Cooley-Tukey butterfly
                    GFT a = block[i];
                    GFT b = block[i + stride];
                    GFT m = gf.mul(b, w);
                    block[i] = gf.add(a, m);
                    block[i + stride] = gf.sub(a, m);
                }
            }
        }
    }

    inline void gs_butterfly(const GFT roots[], GFT block[], size_t block_len, size_t end) const {
        for (size_t stride = block_len / 2, exp_f = 0; stride > 0; stride /= 2, exp_f += 1) {
            for (size_t start = 0; start < end; start += stride * 2) {
                {
                    // Gentleman-Sande butterfly
                    GFT a = block[start];
                    GFT b = block[start + stride];
                    block[start] = gf.add(a, b);
                    block[start + stride] = gf.sub(a, b);
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // Gentleman-Sande butterfly
                    GFT w = roots[(i - start) << exp_f];

                    GFT a = block[i];
                    GFT b = block[i + stride];
                    block[i] = gf.add(a, b);
                    block[i + stride] = gf.mul(gf.sub(a, b), w);
                }
            }
        }
    }

    inline void error_locator(const size_t error_pos_rbo[], size_t error_count, GFT locator_poly[]) const {
        std::fill_n(&locator_poly[0], ecc_len, GFT{0});

        size_t last_error = error_count - 1;

        for (size_t j = 0; j < error_count; ++j) {
            GFT x = GFT{0} + gf.neg(gf.pow(root, error_pos_rbo[j]));
            for (size_t i = last_error - j; i < last_error; ++i) {
                locator_poly[i] = gf.add(locator_poly[i], gf.mul(locator_poly[i + 1], x));
            }
            locator_poly[last_error] = gf.add(locator_poly[last_error], x);
        }
    }

    inline size_t sugiyama(const GFT synd[], GFT a1[], GFT temp5[]) const {
        #define print_vec(v) \
            py::print(#v ":", v##_len, std::vector(&(v)[0], &(v)[block_len]));
        #define print_Vec(v) \
            do { inttr(v); print_vec(v); nttr(v); } while(0)

        #define vec_add(a, b, r) r##_len = _vec_add(a, a##_len, b, b##_len, r)
        #define vec_sub(a, b, r) r##_len = _vec_sub(a, a##_len, b, b##_len, r)
        #define vec_mul(a, b, r) r##_len = _vec_mul(a, a##_len, b, b##_len, r)
        #define vec_div(a, b, r) r##_len = _vec_div(a, a##_len, b, b##_len, r)
        #define vec_copy(a, b) do { \
            std::copy_n(a, ecc_len, b); \
            b##_len = a##_len; \
        } while(0) \

        size_t q_len = 0;
        auto q = &temp5[block_len * 0];

        size_t t_len = 0;
        auto t = &temp5[block_len * 1];

        size_t a1_len;
        std::fill_n(a1, block_len, GFT{0});
        std::fill_n(&temp5[0], block_len * 5, GFT{0});

        auto r1 = &temp5[block_len * 3];
        size_t r1_len = ecc_len - 1;
        std::copy_n(&synd[0], ecc_len - 1, &r1[0]);
        // r1_len = _norm_size(r1, r1_len);
        nttr(r1);

        {
            // Q = R2 // R1
            a1[1] = gf.inv(synd[ecc_len - 1]);
            a1[0] = gf.neg(gf.mul(synd[ecc_len - 2], gf.mul(a1[1], a1[1])));
            a1_len = 2;

            // A1 = A2 - Q * A1
            // a2 = 0
            // a1 = 1
            _vec_neg(a1, q_len, a1);
            nttr(a1);

            // R1 = R2 - Q * R1
            vec_mul(a1, r1, r1);

            inttr(r1);
            r1[ecc_len - 1] = 0;
            r1_len = ecc_len - 1;
            r1_len = _norm_size(r1, r1_len);
            nttr(r1);

            inttr(a1);
            a1_len = _norm_size(a1, a1_len);
            nttr(a1);
        }

        size_t a2_len = 1;
        auto a2 = &temp5[block_len * 4];
        // a2[0] = 1;
        // nttr(a2);
        std::fill_n(&a2[0], block_len, GFT{1});

        auto r2 = &temp5[block_len * 2];
        size_t r2_len = ecc_len;
        std::copy_n(&synd[0], ecc_len, &r2[0]);
        nttr(r2);

        for (size_t i = 1; i < ecc_len / 2 && r1_len > ecc_len / 2; ++i) {
            // Q = R2 // R1

            vec_div(r2, r1, q);

            // t = A2 - Q * A1

            vec_mul(q, a1, t);
            vec_sub(a2, t, t);

            // A2 = A1
            vec_copy(a1, a2);


            // A1 = t if R1.deg() >= len(synds) // 2 else A1
            // TODO: adapt to work with simd
            // if (r1_len > ecc_len)
                vec_copy(t, a1);

            // t = R2 - Q * R1
            vec_mul(q, r1, t);
            vec_sub(r2, t, t);

            // R2 = R1
            vec_copy(r1, r2);

            // R1 = t if R1.deg() >= len(synds) // 2 else R1
            // if (r1_len > ecc_len)
                vec_copy(t, r1);

            inttr(r1);
            r1_len = _norm_size(r1, r1_len);
            nttr(r1);

            inttr(a1);
            a1_len = _norm_size(a1, a1_len);
            nttr(a1);

        }

        // locator = ref.P(GF, [a // GF(A1.x[0]) for a in A1.x])
        inttr(a1);
        GFT a1_0_inv = gf.inv(a1[0]);
        for (size_t i = 0; i < a1_len; ++i)
            a1[i] = gf.mul(a1[i], a1_0_inv);

        std::reverse(&a1[0], &a1[a1_len]);
        return a1_len;
    }

    inline size_t _vec_add(const GFT a[], size_t a_len, const GFT b[], size_t b_len, GFT r[]) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.add(a[i], b[i]);
        return std::max(a_len, b_len);
    }

    inline size_t _vec_sub(const GFT a[], size_t a_len, const GFT b[], size_t b_len, GFT r[]) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.sub(a[i], b[i]);
        return std::max(a_len, b_len);
    }

    inline size_t _vec_neg(const GFT a[], size_t a_len, GFT r[]) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.neg(a[i]);
        return a_len;
    }

    inline size_t _vec_mul(const GFT a[], size_t a_len, const GFT b[], size_t b_len, GFT r[]) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.mul(a[i], b[i]);
        return a_len + b_len - 1;
    }

    inline size_t _norm_size(const GFT r[], size_t r_len) const {
        while (r_len > 0 && r[r_len - 1] == GFT{0})
            --r_len;
        return r_len;
    }

    inline size_t _vec_inv_mod_xn(GFT h[], size_t /*h_len*/, size_t n, GFT a[]) const {
        size_t l = 1;

        // a = P([h.x[0].inv()])
        std::fill_n(&a[0], ecc_len, GFT{0});
        inttr(h);
        GFT h_0 = h[0];
        nttr(h);
        a[0] = gf.inv(h_0);
        nttr(a);

        while (l < n) {
            // A = [a * (GF(2) - a * h) for a, h in zip(A, H)]
            for (size_t i = 0; i < ecc_len; ++i) {
                GFT a_i = a[i];
                GFT ah = gf.mul(a_i, h[i]);
                a[i] = gf.mul(a_i, gf.sub(2, ah));
            }

            // A = A % x**n
            inttr(a);
            std::fill_n(&a[n], ecc_len - n, GFT{0});
            nttr(a);

            l *= 2;
        }

        return n;
    }

    inline size_t _vec_div(GFT f[], size_t f_len, GFT g[], size_t g_len, GFT q[]) const {
        if (f_len < g_len) {
            std::fill_n(q, ecc_len, GFT{0});
            return 0;
        }

        size_t q_len = f_len - g_len + 1;

        // rev_f = P(f.x[::-1])
        // rev_g = P(g.x[::-1])
        inttr(f);
        std::reverse(&f[0], &f[f_len]);
        GFT f_n = f[f_len - 1];
        if (f_len > q_len)
            f[f_len - 1] = 0;
        nttr(f);

        inttr(g);
        std::reverse(&g[0], &g[g_len]);
        nttr(g);

        // rev_g_i = inverse_modulo(rev_g, mod.deg())
        _vec_inv_mod_xn(g, g_len, q_len, q);

        // rev_q = rev_f * rev_g_i % mod
        _vec_mul(f, f_len, q, q_len, q);
        inttr(q);
        std::fill_n(&q[q_len], ecc_len - q_len, GFT{0});
        std::reverse(&q[0], &q[q_len]);
        nttr(q);

        // q = P(rev_q.x[::-1]) * x ** (mod.deg() - rev_q.deg() - 1)
        // std::rotate(&q[0], &q[q_len - 1], &q[q_len]);

        inttr(f);
        f[f_len - 1] = f_n;
        std::reverse(&f[0], &f[f_len]);
        nttr(f);

        inttr(g);
        std::reverse(&g[0], &g[g_len]);
        nttr(g);

        return q_len;
    }

    inline void nttr(GFT block[]) const {
        gs_butterfly(&_roots_v_ecc[0], &block[0], ecc_len, ecc_len);
    }

    inline void inttr(GFT block[]) const {
        ct_butterfly(&_roots_iv_ecc[0], &block[0], ecc_len);
        for (size_t j = 0; j < ecc_len; ++j)
            block[j] = gf.mul(block[j], ecc_len_iv);
    }
};


#undef vec_add
#undef vec_sub
#undef vec_mul
#undef vec_div
#undef vec_copy
