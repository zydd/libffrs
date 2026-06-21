/**************************************************************************
 * ntt.hpp
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
#include <cstdlib>

#include "pygfi16.hpp"
#include "pylogging.hpp"
#include "simd.hpp"


struct NTT {
    ::GFT root;
    size_t ntt_len;

    GFi16 const& gf;
    uint16_t _rbo_shift;
    const size_t block_len;
    const size_t ecc_len;
    const size_t ecc_len_mask;
    const size_t pntt_blocks;
    ::GFT ntt_len_i;
    ::GFT ecc_len_i;
    ::GFT ecc2_len_i;
    std::vector<::GFT> _roots_ntt;
    std::vector<::GFT> _roots_i_ntt;
    std::vector<::GFT> _roots_ecc;
    std::vector<::GFT> _roots_i_ecc;
    std::vector<::GFT> _roots_ecc2;
    std::vector<::GFT> _roots_i_ecc2;
    std::vector<::GFT> _pntt_shift;
    std::vector<::GFT> _pintt_shift;
    std::vector<::GFT> _rbo_ecc;

    inline NTT(GFi16 const& gf, size_t block_len, size_t ecc_len):
        gf(gf),
        block_len(block_len),
        ecc_len(ecc_len),
        ecc_len_mask(ecc_len - 1),
        pntt_blocks(block_len / ecc_len)
    {
        ntt_len = 1 << ffrs::detail::ilog2_ceil(block_len);
        _rbo_shift = 16 - __builtin_ctzl(ntt_len);
        ntt_len_i = gf.inv(ntt_len);
        ecc_len_i = gf.inv(ecc_len);
        ecc2_len_i = gf.inv(ecc_len * 2);

        // root = gf.exp(gf.div(gf.log(1), ntt_len));
        root = gf.pow(gf.primitive, 65536 / ntt_len);
        py_assert(gf.pow(root, ntt_len) == 1, "ntt root of unity sanity check failed");
        py_assert(gf.pow(root, ntt_len / 2) == 65536, "ntt root of unity sanity check failed");

        ::GFT root_i = gf.inv(root);
        _roots_ntt.resize(ntt_len);
        _roots_i_ntt.resize(ntt_len);
        for (size_t i = 0; i < ntt_len; ++i) {
            _roots_ntt[i] = gf.pow(root, i);
            _roots_i_ntt[i] = gf.pow(root_i, i);
        }

        _rbo_ecc.resize(ecc_len);
        for (size_t i = 0; i < ecc_len; ++i) {
            _rbo_ecc[i] = rbo(i * (ntt_len / ecc_len));
        }

        // ::GFT ecc_root = gf.exp(gf.div(gf.log(1), ecc_len));
        // ::GFT ecc_root = gf.pow(root, ntt_len / ecc_len);
        ::GFT ecc_root = gf.pow(gf.primitive, 65536 / ecc_len);
        py_assert(gf.pow(ecc_root, ecc_len) == 1, "ntt root of unity sanity check failed");
        py_assert(gf.pow(ecc_root, ecc_len / 2) == 65536, "ntt root of unity sanity check failed");
        ::GFT ecc_root_i = gf.inv(ecc_root);
        _roots_ecc.resize(ecc_len);
        _roots_i_ecc.resize(ecc_len);
        for (size_t i = 0; i < ecc_len; ++i) {
            _roots_ecc[i] = gf.pow(ecc_root, i);
            _roots_i_ecc[i] = gf.pow(ecc_root_i, i);
        }

        ::GFT ecc2_root = gf.pow(gf.primitive, 65536 / (ecc_len * 2));
        py_assert(gf.pow(ecc2_root, ecc_len * 2) == 1, "ntt root of unity sanity check failed");
        py_assert(gf.pow(ecc2_root, ecc_len) == 65536, "ntt root of unity sanity check failed");
        ::GFT ecc2_root_i = gf.inv(ecc2_root);
        _roots_ecc2.resize(ecc_len * 2);
        _roots_i_ecc2.resize(ecc_len * 2);
        for (size_t i = 0; i < ecc_len * 2; ++i) {
            _roots_ecc2[i] = gf.pow(ecc2_root, i);
            _roots_i_ecc2[i] = gf.pow(ecc2_root_i, i);
        }

        _pntt_shift.resize(block_len);
        for (size_t i = 0; i < block_len; ++i) {
            auto blk = (i / ecc_len) * ecc_len;
            auto j = i - blk;
            _pntt_shift[i] = gf.pow(root, j * rbo(blk));
        }

        _pintt_shift.resize(ntt_len);
        for (size_t i = 0; i < ntt_len; ++i) {
            auto blk = (i / ecc_len) * ecc_len;
            auto j = i - blk;
            _pintt_shift[i] = gf.div(gf.pow(root_i, j * rbo(blk)), ntt_len);
        }
    }

    uint16_t rbo(uint16_t b) const {
        return ffrs::detail::rbo16(b) >> _rbo_shift;
    }

    /**
     * partial NTT
     * computes the first ecc_len symbols of the NTT
     */
    template<typename GFT>
    inline void pntt(GFT *const block) const {
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

    /**
     * partial NTT
     * computes the first `ecc_len` symbols of the NTT
     * this version assumes the last `ecc_len` symbols of `block` are 0
     */
    template<typename GFT>
    inline void pntt_message(GFT *const block) const {
        {
            ct_butterfly(&_roots_ecc[0], &block[0], ecc_len);

            for (size_t j = 0; j < ecc_len; ++j)
                block[j] = gf.mul(block[j], _pntt_shift[j]);
        }

        for (size_t i = 1; i < pntt_blocks - 1; ++i) {
            auto p = &block[i * ecc_len];
            ct_butterfly(&_roots_ecc[0], &p[0], ecc_len);

            for (size_t j = 0; j < ecc_len; ++j)
                p[j] = gf.mul(p[j], _pntt_shift[i * ecc_len + j]);

            for (size_t j = 0; j < ecc_len; ++j)
                block[j] = gf.add(block[j], block[i * ecc_len + j]);
        }
    }

    template<typename GFT>
    inline void pntt_message_residue(GFT *const block) const {
        {
            ct_butterfly_residue(&_roots_ecc[0], &block[0], ecc_len);

            for (size_t j = 0; j < ecc_len; ++j)
                block[j] = gf.mul_residue(gf.mod_p(block[j]), _pntt_shift[j]);
        }

        for (size_t i = 1; i < pntt_blocks - 1; ++i) {
            auto p = &block[i * ecc_len];
            ct_butterfly_residue(&_roots_ecc[0], &p[0], ecc_len);

            for (size_t j = 0; j < ecc_len; ++j)
                p[j] = gf.mul_residue(gf.mod_p(p[j]), _pntt_shift[i * ecc_len + j]);

            for (size_t j = 0; j < ecc_len; ++j)
                block[j] = gf.add_residue(block[j], block[i * ecc_len + j]);
        }
    }

    /**
     * partial iNTT
     * computes the full iNTT from the first `ecc_len` symbols of `block`
     * normal order input, RBO output
     */
    template<typename GFT>
    inline void pintt_ecc(GFT *const ntt_block) const {
        for (size_t i = 1; i < block_len / ecc_len; ++i) {
            std::copy_n(&ntt_block[0], ecc_len, &ntt_block[i * ecc_len]);
            for (size_t j = 0; j < ecc_len; ++j)
                ntt_block[i * ecc_len + j] = gf.mul(ntt_block[i * ecc_len + j], _pintt_shift[i * ecc_len + j]);

            gs_butterfly(&_roots_i_ecc[0], &ntt_block[i * ecc_len], ecc_len, ecc_len);
        }

        for (size_t j = 0; j < ecc_len; ++j)
            ntt_block[j] = gf.mul(ntt_block[j], _pintt_shift[j]);

        gs_butterfly(&_roots_i_ecc[0], &ntt_block[0], ecc_len, ecc_len);
    }

    /**
     * RBO input, normal order output
     */
    template<typename GFT>
    inline void ct_butterfly(const ::GFT *const roots, GFT *const block, size_t ntt_len) const {
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

    template<typename GFT>
    inline void ct_butterfly_residue(const ::GFT *const roots, GFT *const block, size_t ntt_len) const {
        for (size_t stride = 1, exp_f = ntt_len >> 1; stride < ntt_len; stride *= 2, exp_f >>= 1) {
            for (size_t start = 0; start < ntt_len /*input_size*/; start += stride * 2) {
                {
                    // Cooley-Tukey butterfly
                    auto a = block[start];
                    auto b = block[start + stride];
                    block[start] = gf.add_residue(a, b);
                    block[start + stride] = gf.sub_residue(a, b);
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // j = i - start
                    auto w = roots[exp_f * (i - start)];

                    // Cooley-Tukey butterfly
                    auto a = block[i];
                    auto b = block[i + stride];
                    auto m = gf.mul_residue(gf.mod_p(b), w);
                    block[i] = gf.add_residue(a, m);
                    block[i + stride] = gf.sub_residue(a, m);
                }
            }
        }
    }

    /**
     * normal order input, RBO output
     */
    template<typename GFT>
    inline void gs_butterfly(const ::GFT *const roots, GFT *const block, size_t ntt_len, size_t end) const {
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

    template<typename GFT>
    inline void gs_butterfly_residue(const ::GFT *const roots, GFT *const block, size_t ntt_len, size_t end) const {
        for (size_t stride = ntt_len / 2, exp_f = 0; stride > 0; stride /= 2, exp_f += 1) {
            for (size_t start = 0; start < end; start += stride * 2) {
                {
                    // Gentleman-Sande butterfly
                    auto a = block[start];
                    auto b = block[start + stride];
                    block[start] = gf.add_residue(a, b);
                    block[start + stride] = gf.sub_residue(a, b);
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // Gentleman-Sande butterfly
                    auto w = roots[(i - start) << exp_f];

                    auto a = block[i];
                    auto b = block[i + stride];
                    block[i] = gf.add_residue(a, b);
                    block[i + stride] = gf.mul_residue(gf.mod_p(gf.sub_residue(a, b)), w);
                }
            }
        }
    }

    template<typename GFT>
    inline void _reverse_ntt(GFT *const vec, size_t shift) const {
        // check_le_ecc(shift);

        for (size_t i = 1; i < ecc_len; ++i) {
            // rbo(ecc_len - rbo(i))
            auto j = _rbo_ecc[i];
            auto k = _rbo_ecc[ecc_len - j];
            if (i < k)
                std::swap(vec[i], vec[k]);
            vec[i] = gf.mul(vec[i], _roots_i_ecc[j]);
            vec[i] = gf.mul(vec[i], _roots_ecc[shift * j & ecc_len_mask]);
        }
    }

    template<typename GFT>
    inline void _reverse_ntt_vec(GFT *const vec, GFT shift) const {
        for (size_t i = 1; i < ecc_len; ++i) {
            // rbo(ecc_len - rbo(i))
            auto j = _rbo_ecc[i];
            auto k = _rbo_ecc[ecc_len - j];
            if (i < k)
                std::swap(vec[i], vec[k]);
            vec[i] = gf.mul(vec[i], _roots_i_ecc[j]);
            vec[i] = gf.mul(vec[i], gf.gather(&_roots_ecc[0], shift * j & ::GFT(ecc_len_mask)));
        }
    }

    template<typename GFT>
    inline GFT poly_mul_ntt(const GFT *const a, GFT a_len, const GFT *const b, GFT b_len, GFT *const r) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.mul(a[i], b[i]);
        return a_len + b_len - 1;
    }

    template<typename GFT>
    inline GFT poly_inv_mod_xn_rev(GFT *g, GFT g_len, GFT g0, GFT a_len, GFT *const a) const {
        pd_print("poly_inv_mod_xn_rev >", py::arg("sep") = "");
        pd_print_vec("n", &a_len, 1);
        pd_print_vec("rev_g[0]", &g0, 1);
        // check_le_ecc(a_len);

        // a = P([g.x[0].inv()])
        auto a0 = gf.inv(g0);
        std::fill_n(&a[0], ecc_len, a0);

        size_t l = 2;
        auto a_max = vec::max(a_len);
        while (l < a_max) {
            pd_print("inverse modulo", l);
            pd_print_vec_intt("g", g, ecc_len);
            pd_print_vec_intt("a", a, ecc_len);

            // A = [a * (GF(2) - a * g) for a, g in zip(A, G)]
            py_assert(
                vec::max(g_len) + l - 2 <= ecc_len,
                std::to_string(vec::max(g_len)) + " + " + std::to_string(l) + " - 2 > " + std::to_string(ecc_len)
            );
            for (size_t i = 0; i < ecc_len; ++i) {
                GFT a_i = a[i];
                GFT ag = gf.mul(a_i, g[i]);
                a[i] = gf.mul(a_i, gf.sub(GFT{} + 2, ag));
            }

            inttr(a);

            // A = A % x**a_len
            std::fill_n(&a[l], ::GFT(ecc_len - l), GFT{0});

            pd_print_vec("a", a, ecc_len);

            nttr(a);

            l *= 2;
        }

        // last iteration
        {
            pd_print("inverse modulo", l);
            pd_print_vec_intt("g", g, ecc_len);
            pd_print_vec_intt("a", a, ecc_len);

            // // A = [a * (GF(2) - a * g) for a, g in zip(A, G)]
            // for (size_t i = 0; i < ecc_len; ++i) {
            //     GFT a_i = a[i];
            //     GFT ag = gf.mul(a_i, g[i]);
            //     a[i] = gf.mul(a_i, gf.sub(GFT{} + 2, ag));
            // }
            // // Fix aliasing
            // a[0] = a0;
            for (size_t i = 0; i < ecc_len; ++i) {
                GFT ag = gf.mul(a[i], g[i]);
                g[i] = gf.sub(GFT{} + 2, ag);
            }
            pd_print_vec_intt("2 - ag", g, ecc_len);
            for (size_t i = 0; i < ecc_len; ++i)
                g[i] = gf.mul(a[i], g[i]);

            inttr(a);
            inttr(g);
            std::copy_n(&g[l / 2], ecc_len - l / 2, &a[l / 2]);

            // A = A % x**a_len
            std::fill_n(&a[l], ::GFT(ecc_len - l), GFT{0});

            pd_print_vec("a", a, ecc_len);

            // poly_div expects the result to be reversed
            std::reverse(&a[0], &a[a_max]);
            vec::copy_n(&a[0], a_max - a_len, a_len, &a[0], GFT{0});
            vec::fill_n(&a[0], a_len, ::GFT(ecc_len) - a_len, GFT{0});

            pd_print_vec("a_rev", a, ecc_len);

            nttr(a);
        }

        return a_len;
    }

    template<typename GFT>
    inline GFT poly_div_ntt(const GFT *const f, GFT f_len, const GFT *const g, GFT g_len, GFT f0, GFT g0, GFT *const t, GFT *const q) const {
        pd_print("\npoly_div >");
        pd_print_vec("f_len (init)", &f_len, 1);
        pd_print_vec("g_len (init)", &g_len, 1);
        GFT q_len = f_len - g_len + 1;
        pd_print_vec("q_len (init)", &q_len, 1);
        pd_print_vec("f[0]", &f0, 1);
        pd_print_vec("g[-1]", &g0, 1);

        if constexpr (std::is_integral_v<GFT>) {
            py_assert(g_len > 0);
            if (f_len < g_len) {
                std::fill_n(&q[0], ecc_len, GFT{0});
                return 0;
            }
        } else {
            py_assert(!vec::any((g_len == 0) && (f_len > 0)));
            // guard against empty lanes
            q_len &= (g_len != 0) & (f_len >= g_len);
        }

        pd_print_vec("q_len", &q_len, 1);
        // check_le_ecc(q_len);

        std::copy_n(&g[0], ecc_len, &t[0]);
        pd_print_vec_intt("g", t, ecc_len);
        _reverse_ntt_vec(t, g_len);
        pd_print_vec_intt("g_rev", t, ecc_len);

        // auto g_inv_len =
        poly_inv_mod_xn_rev(t, g_len, g0, q_len, q);
        pd_print_vec_intt("g_rev_inv_rev", q, ecc_len);
        pd_print_vec("g_rev_inv_rev_len", &g_inv_len, 1);

        // avoid aliasing: f[0] = 0
        std::copy_n(&f[0], ecc_len, &t[0]);
        pd_print_vec_intt("f", t, ecc_len);
        // if constexpr (std::is_integral_v<GFT>) {
        //     if (f_len > g_len) {
                for (size_t j = 0; j < ecc_len; ++j)
                    t[j] = gf.sub(t[j], f0);
                // inttr(t);
                // vec::fill_n(&t[0], GFT{0}, f_len - q_len + 1, GFT{0});
                // nttr(t);
                pd_print_vec_intt("f (f[0] = 0)", t, ecc_len);
        //     }
        // }

        // t = f * g_rev_inv_rev (shifted) % mod
        poly_mul_ntt(t, f_len, q, q_len, t);

        inttr(t);
        pd_print_vec("f * g_rev_inv_rev", t, ecc_len);
        {

            // if constexpr (std::is_integral_v<GFT>) {
            //     if (f_len > g_len) {
            //         GFT end = (f_len + g_inv_len - 1) & ecc_len_mask;
            //         pd_print_vec("end", &end, 1);
            //         pd_print_vec("q_len", &q_len, 1);

            //         if (end < q_len) {
            //             vec::copy_n(&t[0], GFT{0}, end, &q[0], q_len - end);
            //             pd_print_vec("q1", q, ecc_len);
            //             GFT right = q_len - end;
            //             vec::copy_n(&t[0], ::GFT(ecc_len) - right, right, &q[0], GFT{0});
            //         } else if (end > q_len) {
            //             vec::copy_n(&t[0], end - q_len, q_len, &q[0], GFT{0});
            //         }
            //         q[0] = GFT{0};
            //     } else {
            //         // f_len - q_len = f_len - (f_len - g_len + 1) = g_len - 1
            //         auto right_n = vec::min(q_len, ::GFT(ecc_len) - g_len + 1);
            //         pd_print_vec("right_n", &right_n, 1);
            //         vec::copy_n(&t[0], g_len - 1, right_n, &q[0], GFT{0});
            //         vec::copy_n(&t[0], GFT{0}, q_len - right_n, &q[0], right_n);
            //     }
            // }

            auto right_n = vec::min(q_len, ::GFT(ecc_len) - g_len);
            vec::copy_n(&t[0], g_len, right_n, &q[0], GFT{0});
            vec::copy_n(&t[0], GFT{0}, q_len - right_n, &q[0], right_n);
            vec::fill_n(&q[0], q_len, ::GFT(ecc_len) - q_len, GFT{0});

            vec::fill_n(&q[0], q_len, ::GFT(ecc_len) - q_len, GFT{0});
        }
        pd_print_vec("q", q, ecc_len);

        nttr(q);
        pd_print("poly_div_ntt <\n");

        return q_len;
    }

    /**
     * normal order input, RBO output
     */
    template<typename GFT>
    inline void nttr(GFT *const block) const {
        gs_butterfly(&_roots_ecc[0], &block[0], ecc_len, ecc_len);
    }

    template<typename GFT>
    inline void inttr(GFT *const block) const {
        ct_butterfly(&_roots_i_ecc[0], &block[0], ecc_len);
        for (size_t i = 0; i < ecc_len; ++i)
            block[i] = gf.mul(block[i], ecc_len_i);
    }

    template<typename GFT>
    inline void nttr2(GFT *const block) const {
        gs_butterfly(&_roots_ecc2[0], &block[0], ecc_len * 2, ecc_len * 2);
    }

    template<typename GFT>
    inline void inttr2(GFT *const block) const {
        ct_butterfly(&_roots_i_ecc2[0], &block[0], ecc_len * 2);
        for (size_t i = 0; i < ecc_len * 2; ++i)
            block[i] = gf.mul(block[i], ecc2_len_i);
    }

    /**
     * RBO input, normal order output
     */
    template<typename GFT>
    inline void ntt(GFT *const block) const {
        ct_butterfly(&_roots_ntt[0], &block[0], ntt_len);
    }

    template<typename GFT>
    inline void intt(GFT *const block) const {
        gs_butterfly(&_roots_i_ntt[0], &block[0], ntt_len, ntt_len);
        for (size_t i = 0; i < ntt_len; ++i)
            block[i] = gf.mul(block[i], ntt_len_i);
    }

    template<typename GFT>
    py::bytearray py_ntt(buffer_ro<uint16_t> input) const;
    template<typename GFT>
    py::bytearray py_intt(buffer_ro<uint16_t> input) const;
    template<typename GFT>
    py::bytearray py_nttr(buffer_ro<uint16_t> input) const;
    template<typename GFT>
    py::bytearray py_inttr(buffer_ro<uint16_t> input) const;
    template<typename GFT>
    py::bytearray py_poly_mul(buffer_ro<uint16_t> a, buffer_ro<uint16_t> b) const;
    template<typename GFT>
    py::bytearray py_poly_div(buffer_ro<uint16_t> a, buffer_ro<uint16_t> b) const;
    template<typename GFT>
    py::bytearray py_poly_inv(buffer_ro<uint16_t> a, size_t n) const;

    static void register_class(py::module &m);
};
