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


#include <pybind11/stl.h>

// #define FFRS_CHECK_BOUNDS
// #define FFRS_DEBUG_LONG_DIVISION
// #define FFRS_DEBUG_POLY_DIVISION


template<typename Vec, typename = std::enable_if_t<sizeof(Vec) % sizeof(::GFT) == 0>>
inline void _print_vec(const char *name, const Vec *const v, size_t len) {
        py::print(name, " (", len, ") = [", py::arg("end") = sizeof(Vec)/sizeof(::GFT) > 1 ? "\n" : "", py::arg("sep") = "");
        for (size_t i = 0; i < std::min<size_t>(len, 32); ++i) {
            for (size_t j = 0; j < sizeof(Vec)/sizeof(::GFT); ++j) {
                auto coeff = *(reinterpret_cast<const ::GFT *>(v + i) + j);
                char coeff_fmt[20];
                std::snprintf(coeff_fmt, sizeof(coeff_fmt), sizeof(Vec)/sizeof(::GFT) > 1 ? "%6d" : "%d, ", coeff);
                py::print(coeff_fmt, py::arg("end") = " ");
            }
            if (sizeof(Vec)/sizeof(::GFT) > 1) py::print();
        }
        if (len > 32)
            py::print("   ...");
        py::print(sizeof(Vec)/sizeof(::GFT) > 1 ? "]" : "\b\b\b]");
}


#ifdef FFRS_CHECK_BOUNDS
    inline void _check_size(size_t v, size_t max, const char *location) {
        if (v > max) {
            throw std::runtime_error(
                std::string(location) + ": "
                + std::to_string(v) + " >= " + std::to_string(max)
            );
        }
    }
    template<typename Vec, typename = std::enable_if_t<sizeof(Vec) % sizeof(::GFT) == 0>>
    inline void _check_size(Vec const& v, size_t max, const char *location) {
        auto max_value = *std::max_element(reinterpret_cast<const ::GFT *>(&v), reinterpret_cast<const ::GFT *>(&v + 1));

        if (max_value > max) {
            _print_vec(location, &v, 1);
            throw std::runtime_error(
                std::string(location) + ": "
                + std::to_string(max_value) + " >= " + std::to_string(max)
            );
        }
    }

    #define STRINGIFY(x) #x
    #define TOSTRING(x) STRINGIFY(x)
    #define LOCATION __FILE__ ":" TOSTRING(__LINE__)
    #define check_le_ecc(v) _check_size((v), ecc_len, "check_le_ecc: " LOCATION ": " #v)
#else
    #define check_le_ecc(v)
#endif



#define print_vec(v, len) _print_vec(#v, (v), (len));

#ifdef FFRS_DEBUG_LONG_DIVISION
    #define ld_print_vec(name, v, len) _print_vec(name, v, len)
    #define ld_print_vec_intt(name, v, len) do { inttr(v); _print_vec(name, v, len); nttr(v); } while(0)
    #define ld_print(...) py::print(__VA_ARGS__)
#else
    #define ld_print_vec(...)
    #define ld_print_vec_intt(...)
    #define ld_print(...)
#endif

#ifdef FFRS_DEBUG_POLY_DIVISION
    #define pd_print_vec(name, v, len) _print_vec(name, v, len)
    #define pd_print_vec_intt(name, v, len) do { inttr(v); _print_vec(name, v, len); nttr(v); } while(0)
    #define pd_print(...) py::print(__VA_ARGS__)
#else
    #define pd_print_vec(...)
    #define pd_print_vec_intt(...)
    #define pd_print(...)
#endif


typedef GFT GFTx4 __attribute__((vector_size(4 * sizeof(GFT))));
typedef GFT GFTx8 __attribute__((vector_size(8 * sizeof(GFT))));
typedef GFT GFTx16 __attribute__((vector_size(16 * sizeof(GFT))));


template<int N>
using simd_map_t =
    std::conditional_t<N == 1, GFT,
    std::conditional_t<N == 4, GFTx4,
    std::conditional_t<N == 8, GFTx8,
    std::conditional_t<N == 16, GFTx16,
    void>>>>;

using simd_mask_t = int;

template<typename GFT>
class RSi16vImpl {
public:
    ::GFT root;
    size_t ntt_len;

    inline RSi16vImpl(GFi16 const& gf, size_t block_len, size_t ecc_len):
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

        root = gf.exp(gf.div(gf.log(1), ntt_len));

        if (gf.pow(root, ntt_len) != 1)
            throw std::runtime_error("Root of unity not found for block size");

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

        _ecc_mix.resize(ecc_len);
        _ecc_mix_i.resize(ecc_len);
        auto ecc_mix_w = *reinterpret_cast<::GFT *>(&_roots_ntt[rbo(block_len - ecc_len)]);
        auto ecc_mix_w_i = gf.inv(ecc_mix_w);
        for (size_t i = 0; i < ecc_len; ++i) {
            _ecc_mix[i] = gf.neg(gf.div(gf.pow(ecc_mix_w_i, i), ecc_len));
            _ecc_mix_i[i] = gf.neg(gf.pow(ecc_mix_w, i));
        }

        ::GFT ecc_root = gf.pow(root, ntt_len / ecc_len);
        // ::GFT ecc_root = gf.exp(gf.div(gf.log(1), ecc_len));
        ::GFT ecc_root_i = gf.inv(ecc_root);
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

        _pintt_shift.resize(ntt_len);
        for (size_t i = 0; i < ntt_len; ++i) {
            auto blk = (i / ecc_len) * ecc_len;
            auto j = i - blk;
            _pintt_shift[i] = gf.div(gf.pow(root_i, j * rbo(blk)), ntt_len);
        }
    }

    inline void encode(GFT *const block) const {
        // pntt_message(&block[0]);
        // ecc_mix(&block[0]);
        pntt_message_residue(&block[0]);
        ecc_mix_residue(&block[0]);
    }

    inline void ecc_mix(GFT *const ecc) const {
        for (size_t j = 0; j < ecc_len; ++j)
            ecc[j] = gf.mul(ecc[j], _ecc_mix[j]);

        gs_butterfly(&_roots_i_ecc[0], &ecc[0], ecc_len, ecc_len);
    }

    inline void ecc_mix_residue(GFT *const ecc) const {
        for (size_t j = 0; j < ecc_len; ++j)
            ecc[j] = gf.mul(gf.mod_p(ecc[j]), _ecc_mix[j]);
            // ecc[j] = gf.mul_residue(gf.mod_p(ecc[j]), _ecc_mix[j]);

        gs_butterfly(&_roots_i_ecc[0], &ecc[0], ecc_len, ecc_len);
        // gs_butterfly_residue(&_roots_i_ecc[0], &ecc[0], ecc_len, ecc_len);

        // for (size_t j = 0; j < ecc_len; ++j)
        //     ecc[j] = gf.mod_p(ecc[j]);
    }

    inline void ecc_unmix(GFT *const ecc) const {
        ct_butterfly(&_roots_ecc[0], &ecc[0], ecc_len);

        for (size_t j = 0; j < ecc_len; ++j)
            ecc[j] = gf.mul(ecc[j], _ecc_mix_i[j]);
    }

    /**
     * partial NTT
     * computes the first ecc_len symbols of the NTT
     */
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
     */
    inline void pintt_ecc(GFT *const ntt_block) const {
        for (size_t i = 1; i < ntt_len / ecc_len; ++i) {
            std::copy_n(&ntt_block[0], ecc_len, &ntt_block[i * ecc_len]);
            for (size_t j = 0; j < ecc_len; ++j)
                ntt_block[i * ecc_len + j] = gf.mul(ntt_block[i * ecc_len + j], _pintt_shift[i * ecc_len + j]);

            gs_butterfly(&_roots_i_ecc[0], &ntt_block[i * ecc_len], ecc_len, ecc_len);
        }

        for (size_t j = 0; j < ecc_len; ++j)
            ntt_block[j] = gf.mul(ntt_block[j], _pintt_shift[j]);

        gs_butterfly(&_roots_i_ecc[0], &ntt_block[0], ecc_len, ecc_len);
    }

    inline void repair(GFT *const block, GFT *const temp_ntt1_ecc6) const {
        // temp_ntt1_ecc6 = ntt_len + ecc_len * 6

        // compute synds
        auto synds = &temp_ntt1_ecc6[ecc_len * 0];
        std::copy_n(&block[0], block_len, &temp_ntt1_ecc6[ecc_len * 0]);
        pntt(&synds[ecc_len * 0]);

        // init evaluator_poly with synds
        auto evaluator_poly = &synds[0];

        auto locator_poly = &temp_ntt1_ecc6[ecc_len * 1];
        auto locator_poly_len = sugiyama(&locator_poly[0], &evaluator_poly[0], &temp_ntt1_ecc6[ecc_len * 2]);

        auto evaluator_poly_len = ::GFT(ecc_len) - locator_poly_len + 1;

        auto locator_poly_deriv = &temp_ntt1_ecc6[ecc_len * 2];
        std::copy_n(&locator_poly[0], ecc_len, &locator_poly_deriv[0]);
        _deriv(locator_poly_deriv, ecc_len);
        auto locator_poly_deriv_len = _get_vec_size(&locator_poly_deriv[0], ecc_len);

        auto roots = &temp_ntt1_ecc6[ecc_len * 3];
        find_roots_ntt(
            &block[0],
            &locator_poly[0], locator_poly_len,
            &locator_poly_deriv[0], locator_poly_deriv_len,
            &evaluator_poly[0], evaluator_poly_len,
            &roots[0]
        );
    }

    inline void repair(GFT *const block, const size_t *const error_pos_rbo, size_t error_count, GFT *const temp2) const {
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
    const size_t block_len;
    const size_t ecc_len;
    const size_t ecc_len_mask;
    const size_t pntt_blocks;
    ::GFT ntt_len_i;
    ::GFT ecc_len_i;
    std::vector<::GFT> _roots_ntt;
    std::vector<::GFT> _roots_i_ntt;
    std::vector<::GFT> _roots_ecc;
    std::vector<::GFT> _roots_i_ecc;
    std::vector<::GFT> _ecc_mix;
    std::vector<::GFT> _ecc_mix_i;
    std::vector<::GFT> _pntt_shift;
    std::vector<::GFT> _pintt_shift;
    std::vector<::GFT> _rbo_ecc;

    /**
     * RBO input, normal order output
     */
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
     * Normal order input, RBO output
     */
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

    inline GFT sugiyama(GFT *const a1, GFT *const r1, GFT *const temp_ecc4) const {
        // r1 = synds
        // temp_ecc4 = 4 * ecc_len

        std::fill_n(&temp_ecc4[0], ecc_len * 4, GFT{0});

        GFT a1_len;
        std::fill_n(a1, ecc_len, GFT{0});

        GFT q_len = GFT{0};
        auto q = &temp_ecc4[ecc_len * 0];

        GFT t_len = GFT{0};
        auto t = &temp_ecc4[ecc_len * 1];

        GFT a2_len = GFT{} + 1;
        auto a2 = &temp_ecc4[ecc_len * 2];
        // a2[0] = 1;
        // nttr(a2);
        std::fill_n(&a2[0], ecc_len, GFT{} + 1);

        // TODO test when second-to-last synd is 0
        auto r1_len = _get_vec_size(&r1[0], ecc_len);

        ld_print("\n-----\nInitialize");

        auto r2 = &temp_ecc4[ecc_len * 3];
        std::copy_n(&r1[0], ecc_len, &r2[0]);
        GFT r2_len = r1_len;
        auto r2l = r2[0];

        ld_print_vec("r2_len", &r2_len, 1);
        ld_print_vec("r2/synd", r2, ecc_len);

        {
            // Q = R2 // R1
            // A1 = A2 - Q * A1
            // a2 = 0
            // a1 = 1
            // A1 = -Q
            auto r1h = gather_vec(&r1[0], (r1_len - 1) & ::GFT(ecc_len_mask));
            auto r1h2 = gather_vec(&r1[0], GFT(r1_len - 2) & ::GFT(ecc_len_mask));

            a1[1] = gf.neg(gf.inv(r1h));
            a1[0] = gf.mul(r1h2, gf.mul(a1[1], a1[1]));
            a1_len = GFT{} + 2;

            for (size_t i = 0; i < ecc_len; ++i)
                r1[i] = gf.mul(a1[0], r1[i]);
            for (size_t i = 1; i < ecc_len; ++i)
                r1[i] = gf.add(r1[i], gf.mul(a1[1], r2[i - 1]));

            r1_len = _get_vec_size(r1, ecc_len - 1);
            fill_n_vec(&r1[0], r1_len, ::GFT(ecc_len) - r1_len, GFT{0});

            ld_print_vec("r1_len", &r1_len, 1);
            ld_print_vec("r1", r1, ecc_len);

            ld_print_vec("a1", a1, 2);
            nttr(a1);
            nttr(r2);
        }

        GFT lanes = ~GFT{0};
        if constexpr (!std::is_integral_v<GFT>) {
            lanes = (r1_len != 0);
            a1_len &= lanes;
            a2_len &= lanes;
            r1_len &= lanes;
            r2_len &= lanes;
        }

        ld_print("------\n");
        if (vec_max(r1_len) >= ecc_len / 2) {
            for (size_t i = 1; i < ecc_len / 2; ++i) {
                ld_print("------\nIteration begin");

                // auto r1h = r1[r1_len - 1];
                // TODO: mask gather
                auto r1h = gather_vec(&r1[0], (r1_len - 1) & ::GFT(ecc_len_mask));
                auto r1l = r1[0];
                nttr(r1);

                ld_print_vec_intt("r2", r2, ecc_len);
                ld_print_vec_intt("r1", r1, ecc_len);
                ld_print_vec_intt("a2", a2, ecc_len);
                ld_print_vec_intt("a1", a1, ecc_len);

                // q = r2 // r1
                q_len = _vec_div(r2, r2_len, r1, r1_len, r2l, r1h, t, q);

                ld_print("------");
                ld_print_vec("q_len", &q_len, 1);
                ld_print_vec_intt("q", q, ecc_len);

                // t = q * a1
                t_len = _vec_mul(q, q_len, a1, a1_len, t);
                // q = q * r1
                q_len = _vec_mul(q, q_len, r1, r1_len, q);

                ld_print("------");
                ld_print_vec_intt("r1 * q", q, ecc_len);
                ld_print_vec_intt("a1 * q", t, ecc_len);

                // t = a2 - q * a1
                t_len = _vec_sub(a2, a2_len, t, t_len, t);
                // q = r2 - q * r1
                q_len = _vec_sub(r2, r2_len, q, q_len, q);

                ld_print("------");
                ld_print_vec_intt("r2 - r1 * Q", q, ecc_len);
                ld_print_vec_intt("a2 - a1 * Q", t, ecc_len);
                ld_print();

                // a2 = a1
                std::copy_n(&a1[0], ecc_len, &a2[0]);
                a2_len = a1_len;
                // r2 = r1
                std::copy_n(&r1[0], ecc_len, &r2[0]);
                r2_len = r1_len;
                r2l = r1l;

                {
                    GFT cond = (r1_len >= ::GFT(ecc_len / 2));
                    ld_print_vec("update cond", &cond, 1);

                    if (is_zero(cond)) {
                        inttr(r1);
                        break;
                    }

                    // a1 = t;
                    // r1 = q;
                    copy_n_masked(&t[0], ecc_len, &a1[0], cond);
                    copy_n_masked(&q[0], ecc_len, &r1[0], cond);
                    a1_len = t_len;
                    r1_len = q_len;
                }

                inttr(r1);

                r1_len = _get_vec_size(r1, vec_max(r1_len & lanes));
                ld_print_vec("r1_len", &r1_len, 1);
                ld_print("------\nIteration end");
            }
        }

        inttr(a1);

        a1_len = _get_vec_size(a1, ecc_len);

        ld_print("------\nFinal values");
        ld_print_vec("a1_len", &a1_len, 1);
        ld_print_vec("a1", a1, ecc_len);
        ld_print_vec("r1_len", &r1_len, 1);
        ld_print_vec("r1", r1, ecc_len);

        // locator = ref.P(GF, [a // GF(A1.x[0]) for a in A1.x])
        GFT a1_0_inv = gf.inv(a1[0]);
        // TODO: ecc_len -> max(r1_len)
        for (size_t i = 0; i < ecc_len; ++i)
            a1[i] = gf.mul(a1[i], a1_0_inv);

        // evaluator = ref.P(GF, [a // GF(A1.x[0]) for a in R1.x])
        // TODO: ecc_len -> max(r1_len)
        for (size_t i = 0; i < ecc_len; ++i)
            r1[i] = gf.mul(r1[i], a1_0_inv);

        ld_print_vec("locator", a1, ecc_len);
        ld_print_vec("evaluator", r1, ecc_len);
        ld_print("------\nDone\n");
        return a1_len & lanes;
    }

    inline GFT forney(
        const GFT *const locator_poly_deriv, size_t locator_poly_deriv_len,
        const GFT *const evaluator_poly, size_t evaluator_poly_len,
        ::GFT const& x
    ) const {
        auto x_inv = GFT{} + gf.inv(x);
        auto numerator = _eval(evaluator_poly, evaluator_poly_len, x_inv);
        auto denominator = _eval(locator_poly_deriv, locator_poly_deriv_len, x_inv);
        auto error = gf.mul(gf.div(numerator, denominator), x);
        return error;
    }

    void add_masked(GFT *const vec, GFT const& i, GFT const& value, GFT const& condition) const;
    void assign_masked(GFT& vec, GFT const& value, GFT const& condition) const;
    void copy_n_masked(const GFT *const src, size_t n, GFT *const dst, GFT const& condition) const;

    void copy_n_vec(const GFT *const src, GFT const& src_offset, GFT n, GFT *const dst, GFT const& dst_offset) const {
        if constexpr (std::is_integral_v<GFT>) {
            std::copy_n(&src[src_offset], n, &dst[dst_offset]);
        } else {
            auto max_n = vec_max(n);
            for (size_t i = 0; i < max_n; i++) {
                for (size_t j = 0; j < sizeof(GFT) / sizeof(::GFT); j++)
                    if (i < n[j])
                        dst[i + dst_offset[j]][j] = src[i + src_offset[j]][j];
            }
        }
    }

    void fill_n_vec(GFT *const dst, GFT const& dst_offset, GFT n, GFT const& value) const {
        if constexpr (std::is_integral_v<GFT>) {
            std::fill_n(&dst[dst_offset], n, value);
        } else {
            auto max_n = vec_max(n);
            for (size_t i = 0; i < max_n; i++) {
                for (size_t j = 0; j < sizeof(GFT) / sizeof(::GFT); j++)
                    if (i < n[j])
                        dst[i + dst_offset[j]][j] = value[j];
            }
        }
    }

    bool vec_any(GFT const& a) const {
        if constexpr (std::is_integral_v<GFT>) {
            return a != 0;
        } else {
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                if (a[i] != 0)
                    return true;
            }
            return false;
        }
    }

    ::GFT vec_min(GFT const& a) const {
        if constexpr (std::is_integral_v<GFT>) {
            return a;
        } else {
            auto min = ~::GFT{0};
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                if (a[i] < min)
                    min = a[i];
            }
            return min;
        }
    }

    ::GFT vec_max(GFT const& a) const {
        // *std::max_element(reinterpret_cast<::GFT *>(&root_count), reinterpret_cast<::GFT *>(&root_count + 1));

        if constexpr (std::is_integral_v<GFT>) {
            return a;
        } else {
            ::GFT max = 0;
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                if (a[i] > max)
                    max = a[i];
            }
            return max;
        }
    }

    GFT vec_min(GFT const& a, GFT const& b) const {
        if constexpr (std::is_integral_v<GFT>) {
            return std::min(a, b);
        } else {
            GFT min = GFT{0};
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                min[i] = std::min(a[i], b[i]);
            }
            return min;
        }
    }

    GFT vec_max(GFT const& a, GFT const& b) const {
        if constexpr (std::is_integral_v<GFT>) {
            return std::max(a, b);
        } else {
            GFT max = GFT{0};
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                max[i] = std::max(a[i], b[i]);
            }
            return max;
        }
    }

    void vec_rotate(GFT *const vec, GFT const& n) const {
        if constexpr (std::is_integral_v<GFT>) {
            std::rotate(&vec[0], &vec[n], &vec[ecc_len]);
        } else {
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                auto shift = n[i];
                for (size_t j = 0; j < ecc_len; ++j) {
                    auto coeff = vec[j][i];
                    vec[j][i] = vec[(j + shift) & ecc_len_mask][i];
                    vec[(j + shift) & ecc_len_mask][i] = coeff;
                }
            }
        }
    }

    GFT gather_vec(const GFT *const src, GFT const& i) const {
        if constexpr (std::is_integral_v<GFT>) {
            return src[i];
        } else {
            GFT res = GFT{0};
            for (size_t j = 0; j < sizeof(GFT) / sizeof(::GFT); ++j) {
                check_le_ecc(i[j]);
                res[j] = src[i[j]][j];
            }
            return res;
        }
    }

    inline void error_locator(const size_t *const error_pos_rbo, size_t error_count, GFT *const locator_poly) const {
        std::fill_n(&locator_poly[0], ecc_len, GFT{0});

        size_t last_error = error_count - 1;

        for (size_t j = 0; j < error_count; ++j) {
            GFT x = GFT{} + gf.neg(gf.pow(root, error_pos_rbo[j]));
            for (size_t i = last_error - j; i < last_error; ++i) {
                locator_poly[i] = gf.add(locator_poly[i], gf.mul(locator_poly[i + 1], x));
            }
            locator_poly[last_error] = gf.add(locator_poly[last_error], x);
        }
    }

    inline void find_roots_ntt(
        GFT *const block,
        const GFT *const locator_poly, GFT const& locator_poly_len,
        const GFT *const locator_poly_deriv, GFT const& locator_poly_deriv_len,
        const GFT *const evaluator_poly, GFT const& evaluator_poly_len,
        GFT *const roots
    ) const {
        std::fill_n(&roots[0], ntt_len, GFT{0});
        std::copy_n(&locator_poly[0], ecc_len, &roots[0]);

        // intt(&roots[0]);
        pintt_ecc(&roots[0]);

        GFT lanes = ~GFT{0};
        if constexpr (!std::is_integral_v<GFT>)
            lanes = (locator_poly_len != 0);

        auto locator_poly_deriv_len_max = vec_max(locator_poly_deriv_len & lanes);
        auto evaluator_poly_len_max = vec_max(evaluator_poly_len & lanes);

        for (::GFT i = 0; i < block_len; ++i) {
            auto i_rbo = rbo(i);
            // auto x = gf.pow(root, i_rbo);
            ::GFT x = _roots_ntt[i_rbo];
            GFT root_locations = ((roots[i] | !lanes) == 0);

            if (vec_any(root_locations)) {
                GFT error = forney(
                    &locator_poly_deriv[0], locator_poly_deriv_len_max,
                    &evaluator_poly[0], evaluator_poly_len_max,
                    x
                );

                if constexpr (!std::is_integral_v<GFT>)
                    error &= root_locations;

                // py::print("\nError location:", i);
                // auto error_neg = gf.neg(error);
                // print_vec(&error_neg, 1);
                // auto blk = block[i];
                // print_vec(&blk, 1);

                block[i] = gf.add(block[i], error);
            }
        }
    }

    inline void repair_ntt(GFT *const block_ntt, const GFT *const locator_poly, size_t error_count, GFT *const temp) const {
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

    inline size_t _vec_shift(const GFT *const a, size_t a_len, size_t shift, GFT *const r) const {
        // ::GFT ecc_root = gf.pow(root, ntt_len / ecc_len);

        for (size_t i = 0; i < ecc_len; ++i)
            // r[i] = gf.mul(a[i], gf.pow(ecc_root, shift * _rbo_ecc[i] & ecc_len_mask));
            r[i] = gf.mul(a[i], _roots_ecc[shift * _rbo_ecc[i] & ecc_len_mask]);

        return a_len + shift;
    }

    inline GFT _vec_sub(const GFT *const a, GFT a_len, GFT *const b, GFT b_len, GFT *const r) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.sub(a[i], b[i]);

        return vec_max(a_len, b_len);
    }

    inline size_t _vec_mul(const GFT *const a, size_t a_len, const GFT *const b, size_t b_len, GFT *const r) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.mul(a[i], b[i]);
        return a_len + b_len - 1;
    }

    inline GFT _vec_mul(const GFT *const a, GFT a_len, const GFT *const b, GFT b_len, GFT *const r) const {
        for (size_t i = 0; i < ecc_len; ++i)
            r[i] = gf.mul(a[i], b[i]);
        return a_len + b_len - 1;
    }

    inline size_t _norm_size_rev(GFT *const r, size_t r_len) const {
        check_le_ecc(r_len);
        size_t start = r_len;
        while (start < ecc_len && is_zero(r[start]))
            ++start;

        if (start == ecc_len) {
            start = 0;
            while (start < r_len && is_zero(r[start]))
                ++start;
        }

        std::rotate(&r[0], &r[start], &r[ecc_len]);
        if (r_len == ecc_len && start == 0)
            return r_len;
        else
            return (r_len - start) & ecc_len_mask;
    }

    static simd_mask_t is_zero(GFT const& vec);

    inline GFT _get_vec_size(const GFT *const r, size_t r_len_max) const {
        check_le_ecc(r_len_max);

        GFT len = GFT{} + ::GFT(r_len_max);

        if constexpr (std::is_integral_v<GFT>) {
            while (r_len_max > 0 && is_zero(r[r_len_max - 1]))
                --r_len_max;
            return r_len_max;
        } else {
            GFT searching = (len > 0);
            while (r_len_max > 0) {
                searching &= (r[r_len_max - 1] == 0);
                len += searching;
                searching &= (len > 0);
                --r_len_max;
            }
        }

        check_le_ecc(len);
        return len;
    }

    inline size_t _deriv(GFT *const r, size_t r_len) const {
        check_le_ecc(r_len);

        for (::GFT i = 0; i < r_len - 1; ++i)
            r[i] = gf.mul(r[i + 1], i + 1);
        r[r_len - 1] = GFT{0};
        return r_len - 1;
    }

    template<typename T>
    inline T _eval(const T *const r, size_t r_len, T x) const {
        check_le_ecc(r_len);

        T sum = T{0};
        for (size_t i = r_len - 1; i < r_len; --i)
            sum = gf.add(gf.mul(sum, x), r[i]);
        return sum;
    }

    inline GFT _vec_inv_mod_xn(GFT *const g, GFT g0, GFT a_len, GFT *const a) const {
        pd_print("_vec_inv_mod_xn >", py::arg("sep") = "");
        pd_print_vec("n", &a_len, 1);
        pd_print_vec("rev_g[0]", &g0, 1);
        check_le_ecc(a_len);

        size_t l = 1;

        // a = P([g.x[0].inv()])
        std::fill_n(&a[0], ecc_len, gf.inv(g0));

        auto a_max = vec_max(a_len);
        while (l < a_max) {
            pd_print("inverse modulo", l);

            // A = [a * (GF(2) - a * g) for a, g in zip(A, G)]
            for (size_t i = 0; i < ecc_len; ++i) {
                GFT a_i = a[i];
                GFT ag = gf.mul(a_i, g[i]);
                a[i] = gf.mul(a_i, gf.sub(GFT{} + 2, ag));
            }

            inttr(a);
            // A = A % x**a_len
            // std::fill_n(&a[a_len], ecc_len - a_len, GFT{0});
            fill_n_vec(&a[0], a_len, ::GFT(ecc_len) - a_len, GFT{0});
            pd_print_vec("a", a, ecc_len);

            // Reverse result in last iteration
            if (l * 2 >= a_max) {
                std::reverse(&a[0], &a[a_max]);
                copy_n_vec(&a[0], a_max - a_len, a_len, &a[0], GFT{0});
                fill_n_vec(&a[0], a_len, ::GFT(ecc_len) - a_len, GFT{0});
            }
            nttr(a);

            l *= 2;
        }

        return a_len;
    }

    inline void _reverse_ntt(GFT *const vec, size_t shift) const {
        check_le_ecc(shift);

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

    inline GFT _vec_div(const GFT *const f, GFT f_len, const GFT *const g, GFT g_len, GFT f0, GFT g0, GFT *const t, GFT *const q) const {
        GFT q_len = f_len - g_len + 1;
        pd_print("\n_vec_div >");
        pd_print_vec("q_len (init)", &q_len, 1);

        if constexpr (std::is_integral_v<GFT>) {
            if (f_len <= g_len) {
                std::fill_n(&q[0], ecc_len, GFT{0});
                return 0;
            }
        } else {
            // guard against empty lanes
            q_len &= (g_len != 0) & (f_len > g_len);
        }

        pd_print_vec("q_len", &q_len, 1);
        check_le_ecc(q_len);

        std::copy_n(&g[0], ecc_len, &t[0]);
        pd_print_vec_intt("g", t, ecc_len);
        _reverse_ntt_vec(t, g_len);
        pd_print_vec_intt("g_rev", t, ecc_len);

        _vec_inv_mod_xn(t, g0, q_len, q);
        pd_print_vec_intt("rev_g_inv_rev", q, ecc_len);

        // avoid aliasing: f[0] = 0
        std::copy_n(&f[0], ecc_len, &t[0]);
        pd_print_vec_intt("f", t, ecc_len);
        for (size_t j = 0; j < ecc_len; ++j)
            t[j] = gf.sub(t[j], f0);
        pd_print_vec_intt("f (f[0] = 0)", t, ecc_len);

        // t = f * rev_g_inv_rev (shifted) % mod
        _vec_mul(t, f_len, q, q_len, t);

        inttr(t);
        pd_print_vec("f * rev_g_inv_rev", t, ecc_len);
        {
            auto right_n = vec_min(q_len, ::GFT(ecc_len) - g_len);
            copy_n_vec(&t[0], g_len, right_n, &q[0], GFT{0});
            copy_n_vec(&t[0], GFT{0}, q_len - right_n, &q[0], right_n);
            fill_n_vec(&q[0], q_len, ::GFT(ecc_len) - q_len, GFT{0});
        }
        pd_print_vec("q", q, ecc_len);

        nttr(q);
        pd_print("_vec_div <\n");

        return q_len;
    }

    inline void nttr(GFT *const block) const {
        gs_butterfly(&_roots_ecc[0], &block[0], ecc_len, ecc_len);
    }

    inline void inttr(GFT *const block) const {
        ct_butterfly(&_roots_i_ecc[0], &block[0], ecc_len);
        for (size_t j = 0; j < ecc_len; ++j)
            block[j] = gf.mul(block[j], ecc_len_i);
    }

    inline void ntt(GFT *const block) const {
        ct_butterfly(&_roots_ntt[0], &block[0], ntt_len);
    }

    inline void intt(GFT *const block) const {
        gs_butterfly(&_roots_i_ntt[0], &block[0], ntt_len, ntt_len);
        for (size_t j = 0; j < ntt_len; ++j)
            block[j] = gf.mul(block[j], ntt_len_i);
    }
};


template <size_t W>
class RSi16v<W>::Impl : public RSi16vImpl<simd_map_t<W>> {
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
void RSi16v<W>::encode(GFT *const block) const {
    d->encode(reinterpret_cast<simd_map_t<W> *>(block));
}


template <size_t W>
void RSi16v<W>::ecc_mix(GFT *const block) const {
    d->ecc_mix(reinterpret_cast<simd_map_t<W> *>(block));
}


template <size_t W>
void RSi16v<W>::pntt(GFT *const block) const {
    d->pntt(reinterpret_cast<simd_map_t<W> *>(block));
}


template <size_t W>
uint16_t RSi16v<W>::rbo(uint16_t v) const {
    return d->rbo(v);
}


template <size_t W>
void RSi16v<W>::repair(GFT *const block, const size_t *const error_pos_rbo, size_t error_count, GFT *const temp2) const {
    d->repair(
        reinterpret_cast<simd_map_t<W> *>(block),
        error_pos_rbo,
        error_count,
        reinterpret_cast<simd_map_t<W> *>(temp2)
    );
}


template <size_t W>
void RSi16v<W>::repair(GFT *const block, GFT *const temp1_ecc6) const {
    d->repair(
        reinterpret_cast<simd_map_t<W> *>(block),
        reinterpret_cast<simd_map_t<W> *>(temp1_ecc6)
    );
}
