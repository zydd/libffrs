/**************************************************************************
 * rsi16md_impl.hpp
 *
 * Copyright 2026 Gabriel Machado
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
#include "rsi16md.h"


#define print_vec(v) \
    py::print(#v ":", v##_len, std::vector(&(v)[0], &(v)[v##_len]));
#define print_Vec(v) \
    do { inttr(v); print_vec(v); nttr(v); } while(0)

#define vec_sub(a, b, r) r##_len = _vec_sub(a, a##_len, b, b##_len, r)
#define vec_mul(a, b, r) r##_len = _vec_mul(a, a##_len, b, b##_len, r)
#define vec_div(a, b, r) r##_len = _vec_div(a, a##_len, b, b##_len, r)
#define vec_copy(a, b) \
    do { std::copy_n(a, ecc_len, b); b##_len = a##_len; } while(0)


typedef uint32_t GFTx4 __attribute__((vector_size(4 * sizeof(uint32_t))));
typedef uint32_t GFTx8 __attribute__((vector_size(8 * sizeof(uint32_t))));
typedef uint32_t GFTx16 __attribute__((vector_size(16 * sizeof(uint32_t))));


template<int N>
using simd_map_t =
    std::conditional_t<N == 1, GFT,
    std::conditional_t<N == 4, GFTx4,
    std::conditional_t<N == 8, GFTx8,
    std::conditional_t<N == 16, GFTx16,
    void>>>>;


template<typename GFT>
class RSi16vImpl {
public:
    uint32_t root;
    size_t ntt_len;

    inline RSi16vImpl(GFi16 const& gf, size_t block_len, size_t ecc_len):
        gf(gf),
        block_len(block_len),
        ecc_len(ecc_len),
        pntt_blocks(block_len / ecc_len)
    {
        ntt_len = 1 << ffrs::detail::ilog2_ceil(block_len);
        _rbo_shift = 16 - __builtin_ctzl(ntt_len);
        ntt_len_i = gf.inv(ntt_len);
        ecc_len_i = gf.inv(ecc_len);

        root = gf.exp(gf.div(gf.log(1), ntt_len));

        if (gf.pow(root, ntt_len) != 1)
            throw std::runtime_error("Root of unity not found for block size");

        uint32_t root_i = gf.inv(root);
        _roots_block.resize(block_len);
        _roots_i_block.resize(block_len);
        for (size_t i = 0; i < block_len; ++i) {
            _roots_block[i] = gf.pow(root, i);
            _roots_i_block[i] = gf.pow(root_i, i);
        }

        _ecc_mix.resize(ecc_len);
        _ecc_mix_i.resize(ecc_len);
        auto ecc_mix_w = *reinterpret_cast<uint32_t *>(&_roots_block[rbo(block_len - ecc_len)]);
        auto ecc_mix_w_i = gf.inv(ecc_mix_w);
        for (size_t i = 0; i < ecc_len; ++i) {
            _ecc_mix[i] = gf.neg(gf.div(gf.pow(ecc_mix_w_i, i), ecc_len));
            _ecc_mix_i[i] = gf.neg(gf.pow(ecc_mix_w, i));
        }

        uint32_t ecc_root = gf.pow(root, ntt_len / ecc_len);
        // uint32_t ecc_root = gf.exp(gf.div(gf.log(1), ecc_len));
        uint32_t ecc_root_i = gf.inv(ecc_root);
        _roots_ecc.resize(ecc_len);
        _roots_i_ecc.resize(ecc_len);
        for (size_t i = 0; i < ecc_len; ++i) {
            _roots_ecc[i] = gf.pow(ecc_root, i);
            _roots_i_ecc[i] = gf.pow(ecc_root_i, i);
        }

        _pntt_shift.resize(block_len);
        for (size_t i = 0; i < block_len; ++i) {
            auto blk = (i / ecc_len) * ecc_len;
            auto j = i - blk;
            _pntt_shift[i] = gf.pow(root, j * rbo(blk));
        }
    }

    inline void encode(uint32_t _block[]) const {
        auto block = reinterpret_cast<GFT *>(_block);
        pntt(&block[0]);
        ecc_mix(&block[0]);
    }

    inline void ecc_mix(GFT ecc[]) const {
        for (size_t j = 0; j < ecc_len; ++j)
            ecc[j] = gf.mul(ecc[j], _ecc_mix[j]);

        gs_butterfly(&_roots_i_ecc[0], &ecc[0], ecc_len, ecc_len);
    }

    inline void ecc_unmix(GFT ecc[]) const {
        ct_butterfly(&_roots_ecc[0], &ecc[0], ecc_len);

        for (size_t j = 0; j < ecc_len; ++j)
            ecc[j] = gf.mul(ecc[j], _ecc_mix_i[j]);
    }

    inline void pntt(GFT block[]) const {
        // partial ntt
        // computes the first n symbols of the ntt of size m, given m is divisible by n
        {
            ct_butterfly(&_roots_ecc[0], &block[0], ecc_len);

            for (size_t j = 0; j < ecc_len; ++j)
                block[j] = gf.mul(block[j], _pntt_shift[j]);
        }

        for (size_t i = 1; i < pntt_blocks; ++i) {
            auto p = &block[i * ecc_len];
            ct_butterfly(&_roots_ecc[0], &p[0], ecc_len);

            for (size_t j = 0; j < ecc_len; ++j)
                p[j] = gf.mul(p[j], _pntt_shift[i * ecc_len + j]);

            for (size_t j = 0; j < ecc_len; ++j)
                block[j] = gf.add(block[j], block[i * ecc_len + j]);
        }
    }

    inline void repair(GFT block[], GFT temp1_ecc6[]) const {
        // auto block = reinterpret_cast<GFT *>(_block);
        // auto temp1_ecc6 = reinterpret_cast<GFT *>(_temp1_ecc6);
        // temp1_ecc6 = max(block_len, ecc_len * 6)

        // init evaluator_poly with syndromes
        auto evaluator_poly = &temp1_ecc6[0];
        std::copy_n(&block[0], block_len, &evaluator_poly[0]);
        pntt(evaluator_poly);

        auto locator_poly = &temp1_ecc6[ecc_len];
        auto locator_poly_len = sugiyama(&locator_poly[0], &evaluator_poly[0], &temp1_ecc6[ecc_len * 2]);
        auto evaluator_poly_len = ecc_len - locator_poly_len + 1;

        auto error_locations = &temp1_ecc6[ecc_len * 2];
        auto error_locations_rbo = &temp1_ecc6[ecc_len * 3];
        auto error_count = find_roots(&locator_poly[0], locator_poly_len, &error_locations[0], &error_locations_rbo[0]);

        locator_poly_len = _deriv(locator_poly, locator_poly_len);

        forney(
            &block[0],
            &locator_poly[0], locator_poly_len,
            &evaluator_poly[0], evaluator_poly_len,
            &error_locations[0], &error_locations_rbo[0], error_count
        );
    }

    inline void repair(GFT block[], const size_t error_pos_rbo[], size_t error_count, GFT temp2[]) const {
        // auto block = reinterpret_cast<GFT *>(_block);
        // auto temp2 = reinterpret_cast<GFT *>(_temp2);

        ntt(block);

        auto locator_poly = &temp2[block_len];
        error_locator(error_pos_rbo, error_count, locator_poly);
        repair_ntt(block, locator_poly, error_count, temp2);

        intt(block);
    }

    uint16_t rbo(uint16_t b) const {
        return ffrs::detail::rbo16(b) >> _rbo_shift;
    }

protected:
    GFi16 const& gf;
    uint16_t _rbo_shift;
    size_t block_len;
    size_t ecc_len;
    size_t pntt_blocks;
    uint32_t ntt_len_i;
    uint32_t ecc_len_i;
    std::vector<uint32_t> _roots_block;
    std::vector<uint32_t> _roots_i_block;
    std::vector<uint32_t> _roots_ecc;
    std::vector<uint32_t> _roots_i_ecc;
    std::vector<uint32_t> _ecc_mix;
    std::vector<uint32_t> _ecc_mix_i;
    std::vector<uint32_t> _pntt_shift;

    /**
     * RBO input, normal order output
     */
    inline void ct_butterfly(const uint32_t roots[], GFT block[], size_t ntt_len) const {
        for (size_t stride = 1, exp_f = ntt_len >> 1; stride < ntt_len; stride *= 2, exp_f >>= 1) {
            for (size_t start = 0; start < ntt_len /*input_size*/; start += stride * 2) {
                {
                    // Cooley-Tukey butterfly
                    auto a = block[start];
                    auto b = block[start + stride];
                    block[start] = gf.add(a, b);
                    block[start + stride] = gf.sub(a, b);
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // j = i - start
                    auto w = roots[exp_f * (i - start)];

                    // Cooley-Tukey butterfly
                    auto a = block[i];
                    auto b = block[i + stride];
                    auto m = gf.mul(b, w);
                    block[i] = gf.add(a, m);
                    block[i + stride] = gf.sub(a, m);
                }
            }
        }
    }

    /**
     * Normal order input, RBO output
     */
    inline void gs_butterfly(const uint32_t roots[], GFT block[], size_t ntt_len, size_t end) const {
        for (size_t stride = ntt_len / 2, exp_f = 0; stride > 0; stride /= 2, exp_f += 1) {
            for (size_t start = 0; start < end; start += stride * 2) {
                {
                    // Gentleman-Sande butterfly
                    auto a = block[start];
                    auto b = block[start + stride];
                    block[start] = gf.add(a, b);
                    block[start + stride] = gf.sub(a, b);
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // Gentleman-Sande butterfly
                    auto w = roots[(i - start) << exp_f];

                    auto a = block[i];
                    auto b = block[i + stride];
                    block[i] = gf.add(a, b);
                    block[i + stride] = gf.mul(gf.sub(a, b), w);
                }
            }
        }
    }

    inline size_t sugiyama(GFT a1[], GFT r1[], GFT temp_ecc4[]) const {
        // r1 = synds
        // temp_ecc4 = 4 * ecc_len

        std::fill_n(&temp_ecc4[0], ecc_len * 4, GFT{0});

        size_t a1_len;
        std::fill_n(a1, ecc_len, GFT{0});

        // TODO test when second-to-last synd is 0
        size_t r1_len = _norm_size(r1, ecc_len - 1);

        size_t q_len = 0;
        auto q = &temp_ecc4[ecc_len * 0];

        size_t t_len = 0;
        auto t = &temp_ecc4[ecc_len * 1];

        size_t a2_len = 1;
        auto a2 = &temp_ecc4[ecc_len * 2];
        // a2[0] = 1;
        // nttr(a2);
        std::fill_n(&a2[0], ecc_len, GFT{1});

        auto r2 = &temp_ecc4[ecc_len * 3];
        size_t r2_len = ecc_len;
        std::copy_n(&r1[0], ecc_len, &r2[0]);

        {
            // Q = R2 // R1
            // A1 = A2 - Q * A1
            // a2 = 0
            // a1 = 1
            // A1 = -Q
            a1[1] = gf.neg(gf_inv(r1[ecc_len - 1]));
            a1[0] = gf.mul(r1[ecc_len - 2], gf.mul(a1[1], a1[1]));
            a1_len = 2;

            std::fill_n(&r1[r1_len], ecc_len - r1_len, GFT{0});
            for (size_t i = 0; i < r1_len; ++i)
                r1[i] = gf.mul(a1[0], r1[i]);
            for (size_t i = 1; i < r1_len; ++i)
                r1[i] = gf.add(r1[i], gf.mul(a1[1], r2[i - 1]));

            r1_len = _norm_size(r1, r1_len + 1);
            a1_len = _norm_size(a1, a1_len);

            nttr(a1);
            // nttr(r1); // div ntt elision
            // nttr(r2); // div ntt elision
        }

        for (size_t i = 1; i < ecc_len / 2 && r1_len > ecc_len / 2; ++i) {
            // q = r2 // r1
            vec_div(r2, r1, q);

            // t = q * a1
            vec_mul(q, a1, t);
            // q = q * r1
            vec_mul(q, r1, q);

            // t = a2 - q * a1
            vec_sub(a2, t, t);
            // t = r2 - q * r1
            vec_sub(r2, q, q);

            // a2 = a1
            vec_copy(a1, a2);
            // r2 = r1
            vec_copy(r1, r2);

            // TODO: adapt to work with simd
            // if (r1_len > ecc_len)
                vec_copy(t, a1);
                vec_copy(q, r1);

            inttr(a1);
            inttr(r1);
            inttr(r2); // div ntt elision
            r1_len = _norm_size(r1, r1_len);
            a1_len = _norm_size(a1, a1_len);
            nttr(a1);
            // nttr(r1); // div ntt elision
        }

        // locator = ref.P(GF, [a // GF(A1.x[0]) for a in A1.x])
        inttr(a1);
        GFT a1_0_inv = gf_inv(a1[0]);
        for (size_t i = 0; i < a1_len; ++i)
            a1[i] = gf.mul(a1[i], a1_0_inv);

        // std::reverse(&a1[0], &a1[a1_len]);

        // evaluator = ref.P(GF, [a // GF(A1.x[0]) for a in R1.x])
        // inttr(r1); // div ntt elision
        for (size_t i = 0; i < r1_len; ++i)
            r1[i] = gf.mul(r1[i], a1_0_inv);

        return a1_len;
    }

    inline size_t forney(
        GFT block[],
        const GFT locator_poly_deriv[], size_t locator_poly_deriv_len,
        const GFT evaluator_poly[], size_t evaluator_poly_len,
        const GFT error_pos[], const GFT error_pos_rbo[], size_t error_count
    ) const {
        for (size_t i = 0; i < error_count; ++i) {
            auto x = error_pos_rbo[i];
            auto x_inv = gf_inv(x);
            auto numerator = _eval(evaluator_poly, evaluator_poly_len, x_inv);
            auto denominator = _eval(locator_poly_deriv, locator_poly_deriv_len, x_inv);

            auto error = gf.mul(gf_div(numerator, denominator), x);
            auto is_error = (error_pos[i] < 0x10001);

            // block[error_pos[i]] = gf.add(gather(block, error_pos[i]), error & is_error);
            scatter(&block[0], error_pos[i], gf.add(gather(&block[0], error_pos[i]), gf.mul(error, is_error)));
        }
        return error_count;
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

    inline size_t find_roots(const GFT locator_poly[], size_t locator_poly_len, GFT roots[], GFT roots_rbo[]) const {
        size_t root_count = 0;
        for (uint32_t i = 0; i < block_len; ++i) {
            auto x = gf.pow(root, rbo(i));
            auto x_inv = gf.inv(x);
            auto is_zero = (_eval(locator_poly, locator_poly_len, GFT{x_inv}) == 0);

            if constexpr (std::is_integral_v<GFT>) {
                if (is_zero) {
                    roots[root_count] = i;
                    roots_rbo[root_count] = x;
                    ++root_count;
                }
            } else {
                for (uint32_t j = 0; i < sizeof(GFT) / sizeof(uint32_t); ++i) {
                    if (is_zero[j]) {
                        roots[root_count] = (is_zero & GFT{i}) | !is_zero;
                        roots_rbo[root_count] = (is_zero & GFT{x}) | !is_zero;
                        ++root_count;
                        break;
                    }
                }
            }
        }
        return root_count;
    }

    inline void repair_ntt(GFT block_ntt[], const GFT locator_poly[], size_t error_count, GFT temp[]) const {
        // copy synds
        std::copy_n(&block_ntt[0], ecc_len, &temp[0]);

        for (size_t j = ecc_len; j < ntt_len; ++j) {
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
        //     block[i] = gf.div(block[i], ntt_len);
        // }
    }

    inline GFT gf_inv(GFT const& a) const {
        if constexpr (std::is_integral_v<GFT>) {
            return gf.inv(a);
        } else {
            // return gf._exp[gf.field_elements-1 - _log[a]];
            GFT i = GFT{0x10000} - gather(&gf.log(0), a);
            return gather(&gf.exp(0), i);
        }
    }

    inline GFT gf_div(GFT const& a, GFT const& b) const {
        if constexpr (std::is_integral_v<GFT>) {
            return gf.div(a, b);
        } else {
            auto q = gather(&gf.log(0), a) + 0x10000 - gather(&gf.log(0), b);
            q -= (q >= 0x10000) & 0x10000;
            return gather(&gf.exp(0), q);
        }
    }

    template<typename B>
    inline GFT gf_pow(B const& b, GFT const& e) const {
        if constexpr (std::is_integral_v<GFT>) {
            return gf.pow(b, e);
        } else {
            // return _exp[(_log[b] * e) % (gf.field_elements - 1)];
            GFT p;
            if constexpr (std::is_integral_v<B>) {
                p = GFT{gf.log(b)};
            } else {
                p = gather(&gf.log(0), b);
            }

            p = (p * e) & 0xffff;
            return gather(&gf.exp(0), p);
        }
    }

    inline size_t _vec_sub(const GFT a[], size_t a_len, const GFT b[], size_t b_len, GFT r[]) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.sub(a[i], b[i]);
        return std::max(a_len, b_len);
    }

    inline size_t _vec_mul(const GFT a[], size_t a_len, const GFT b[], size_t b_len, GFT r[]) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.mul(a[i], b[i]);
        return a_len + b_len - 1;
    }

    inline size_t _norm_size(const GFT r[], size_t r_len) const {
        if constexpr (std::is_integral_v<GFT>) {
             while (r_len > 0 && r[r_len - 1] == GFT{0})
                --r_len;
        } else {
            while (r_len > 0) {
                for (size_t j = 0; j < sizeof(GFT) / sizeof(uint32_t); ++j) {
                    if (r[r_len - 1][j] != 0)
                        return r_len;
                }
                --r_len;
            }
        }
        return r_len;
    }

    inline size_t _deriv(GFT r[], size_t r_len) const {
        for (uint32_t i = 0; i < r_len - 1; ++i)
            r[i] = gf.mul(r[i + 1], i + 1);
        r[r_len - 1] = GFT{0};
        return r_len - 1;
    }

    template<typename T>
    inline T _eval(const T r[], size_t r_len, T x) const {
        T sum = T{0};
        for (size_t i = r_len - 1; i < r_len; --i)
            sum = gf.add(gf.mul(sum, x), r[i]);
        return sum;
    }

    inline size_t _vec_inv_mod_xn(GFT h[], GFT h0, size_t n, GFT a[]) const {
        size_t l = 1;

        // a = P([h.x[0].inv()])
        std::fill_n(&a[0], ecc_len, gf_inv(h0));

        while (l < n) {
            // A = [a * (GF(2) - a * h) for a, h in zip(A, H)]
            for (size_t i = 0; i < ecc_len; ++i) {
                GFT a_i = a[i];
                GFT ah = gf.mul(a_i, h[i]);
                a[i] = gf.mul(a_i, gf.sub(GFT{2}, ah));
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
        // inttr(f); // div ntt elision
        std::reverse(&f[0], &f[f_len]);
        GFT f_n = f[f_len - 1];
        if (f_len > q_len)
            f[f_len - 1] = GFT{0};
        nttr(f);

        // inttr(g); // div ntt elision
        std::reverse(&g[0], &g[g_len]);
        auto g0 = g[0];
        nttr(g);

        // rev_g_i = inverse_modulo(rev_g, mod.deg())
        _vec_inv_mod_xn(g, g0, q_len, q);

        // rev_q = rev_f * rev_g_i % mod
        _vec_mul(f, f_len, q, q_len, q);
        inttr(q);
        std::fill_n(&q[q_len], ecc_len - q_len, GFT{0});
        std::reverse(&q[0], &q[q_len]);
        nttr(q);

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
        gs_butterfly(&_roots_ecc[0], &block[0], ecc_len, ecc_len);
    }

    inline void inttr(GFT block[]) const {
        ct_butterfly(&_roots_i_ecc[0], &block[0], ecc_len);
        for (size_t j = 0; j < ecc_len; ++j)
            block[j] = gf.mul(block[j], ecc_len_i);
    }

    inline void ntt(GFT block[]) const {
        ct_butterfly(&_roots_block[0], &block[0], ntt_len);
    }

    inline void intt(GFT block[]) const {
        gs_butterfly(&_roots_i_block[0], &block[0], ntt_len, ntt_len);
        for (size_t j = 0; j < ntt_len; ++j)
            block[j] = gf.mul(block[j], ntt_len_i);
    }

    GFT gather(const uint32_t vec[], GFT const& i) const;
    template<typename T=GFT, typename=std::enable_if_t<!std::is_integral_v<T>>>
    GFT gather(const GFT vec[], GFT const& i) const;

    void scatter(uint32_t vec[], GFT const& i, GFT const& v) const;
    template<typename T=GFT, typename=std::enable_if_t<!std::is_integral_v<T>>>
    void scatter(GFT vec[], GFT const& i, GFT const& v) const;
};


template <size_t W>
struct RSi16v<W>::Impl : RSi16vImpl<simd_map_t<W>> {
    using RSi16vImpl<simd_map_t<W>>::RSi16vImpl;
};


template <size_t W>
RSi16v<W>::RSi16v(GFi16 const& gf, size_t block_size, size_t ecc_len):
    d(new Impl(gf, block_size, ecc_len)),
    root(d->root),
    ntt_len(d->ntt_len)
{ }


template <size_t W>
RSi16v<W>::~RSi16v() {
    delete d;
    d = nullptr;
}


template <size_t W>
void RSi16v<W>::encode(GFT block[]) const {
    d->encode(reinterpret_cast<GFT *>(block));
}


template <size_t W>
void RSi16v<W>::ecc_mix(GFT block[]) const {
    d->ecc_mix(reinterpret_cast<simd_map_t<W> *>(block));
}


template <size_t W>
void RSi16v<W>::pntt(GFT block[]) const {
    d->pntt(reinterpret_cast<simd_map_t<W> *>(block));
}


template <size_t W>
uint16_t RSi16v<W>::rbo(uint16_t v) const {
    return d->rbo(v);
}


template <size_t W>
void RSi16v<W>::repair(GFT block[], const size_t error_pos_rbo[], size_t error_count, GFT temp2[]) const {
    d->repair(
        reinterpret_cast<simd_map_t<W> *>(block),
        error_pos_rbo,
        error_count,
        reinterpret_cast<simd_map_t<W> *>(temp2)
    );
}


template <size_t W>
void RSi16v<W>::repair(GFT block[], GFT temp2[]) const {
    d->repair(
        reinterpret_cast<simd_map_t<W> *>(block),
        reinterpret_cast<simd_map_t<W> *>(temp2)
    );
}


#undef print_vec
#undef print_Vec
#undef vec_sub
#undef vec_mul
#undef vec_div
#undef vec_copy
