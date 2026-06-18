/**************************************************************************
 * rsi16v_impl.hpp
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
#include "rsi16v.hpp"
#include "ntt.hpp"
#include "pylogging.hpp"
#include "simd.hpp"

#include <pybind11/stl.h>

// #define FFRS_CHECK_BOUNDS


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


template<size_t W>
RSi16v<W>::RSi16v(GFi16 const& gf, size_t block_len, size_t ecc_len):
    gf(gf),
    ntt(gf, block_len, ecc_len),
    root(ntt.root),
    ntt_len(ntt.ntt_len),
    block_len(block_len),
    ecc_len(ecc_len),
    ecc_len_mask(ecc_len - 1)
{
    _ecc_mix.resize(ecc_len);
    _ecc_mix_i.resize(ecc_len);
    auto ecc_mix_w = *reinterpret_cast<const ::GFT *>(&ntt._roots_ntt[ntt.rbo(block_len - ecc_len)]);
    auto ecc_mix_w_i = gf.inv(ecc_mix_w);
    for (size_t i = 0; i < ecc_len; ++i) {
        _ecc_mix[i] = gf.neg(gf.div(gf.pow(ecc_mix_w_i, i), ecc_len));
        _ecc_mix_i[i] = gf.neg(gf.pow(ecc_mix_w, i));
    }
}


template<size_t W>
RSi16v<W>::~RSi16v() { }


template<size_t W>
void RSi16v<W>::_encode(GFT *const block) const {
    ntt.pntt_message_residue(&block[0]);
    _mix_ecc_residue(&block[0]);
}


template<size_t W>
void RSi16v<W>::_mix_ecc(GFT *const ecc) const {
    for (size_t j = 0; j < ecc_len; ++j)
        ecc[j] = gf.mul(ecc[j], _ecc_mix[j]);

    ntt.gs_butterfly(&ntt._roots_i_ecc[0], &ecc[0], ecc_len, ecc_len);
}

template<size_t W>
void RSi16v<W>::_mix_ecc_residue(GFT *const ecc) const {
    for (size_t j = 0; j < ecc_len; ++j)
        ecc[j] = gf.mul(gf.mod_p(ecc[j]), _ecc_mix[j]);
        // ecc[j] = gf.mul_residue(gf.mod_p(ecc[j]), _ecc_mix[j]);

    ntt.gs_butterfly(&ntt._roots_i_ecc[0], &ecc[0], ecc_len, ecc_len);
    // ntt.gs_butterfly_residue(&ntt._roots_i_ecc[0], &ecc[0], ecc_len, ecc_len);

    // for (size_t j = 0; j < ecc_len; ++j)
    //     ecc[j] = gf.mod_p(ecc[j]);
}


template<size_t W>
uint16_t RSi16v<W>::rbo(uint16_t v) const {
    return ntt.rbo(v);
}


template<size_t W>
void RSi16v<W>::_repair(GFT *const block, GFT *const temp_ntt1_ecc6) const {
    r_print("repair unknown >");
    // temp_ntt1_ecc6 = ntt_len + ecc_len * 6

    // compute synds
    auto synds = &temp_ntt1_ecc6[ecc_len * 0];
    std::copy_n(&block[0], block_len, &synds[0]);
    ntt.pntt(&synds[0]);
    r_print_vec("synds", synds, ecc_len);

    // stop if all synds are zero
    if (std::all_of(&synds[0], &synds[ecc_len], [](const GFT& v) { return vec<GFT>::max(v) == 0; })) {
        r_print("no errors detected");
        return;
    }

    // init evaluator_poly with synds
    auto evaluator_poly = &synds[0];

    auto locator_poly = &temp_ntt1_ecc6[ecc_len * 1];
    auto locator_poly_len = _sugiyama(&locator_poly[0], &evaluator_poly[0], &temp_ntt1_ecc6[ecc_len * 2]);
    r_print_vec("locator_poly", locator_poly, ecc_len);
    r_print_vec("locator_poly_len", &locator_poly_len, 1);
    r_print_vec("evaluator_poly", evaluator_poly, ecc_len);

    auto evaluator_poly_len = ::GFT(ecc_len) - locator_poly_len + 1;
    if constexpr (!std::is_integral_v<GFT>)
        evaluator_poly_len &= (locator_poly_len != 0);

    auto locator_poly_deriv = &temp_ntt1_ecc6[ecc_len * 2];
    std::copy_n(&locator_poly[0], ecc_len, &locator_poly_deriv[0]);
    _deriv(locator_poly_deriv, ecc_len);
    r_print_vec("locator_poly_deriv", locator_poly_deriv, ecc_len);
    auto locator_poly_deriv_len = _get_vec_size(&locator_poly_deriv[0], ecc_len);

    auto roots = &temp_ntt1_ecc6[ecc_len * 3];
    _find_roots_ntt(
        &block[0],
        &locator_poly[0], locator_poly_len,
        &locator_poly_deriv[0], locator_poly_deriv_len,
        &evaluator_poly[0], evaluator_poly_len,
        &roots[0]
    );
    r_print("repair <");
}


template<size_t W>
void RSi16v<W>::_repair(GFT *const block, const size_t *const error_pos_rbo, size_t error_count, GFT *const temp_ntt1_ecc6) const {
    r_print("repair known >");
    // temp_ntt1_ecc6 = ntt_len + ecc_len * 6

    // compute synds
    auto synds = &temp_ntt1_ecc6[ecc_len * 0];
    std::copy_n(&block[0], block_len, &synds[0]);
    ntt.pntt(&synds[0]);
    r_print_vec("synds", synds, ecc_len);

    // stop if all synds are zero
    if (std::all_of(&synds[0], &synds[ecc_len], [](const GFT& v) { return vec<GFT>::max(v) == 0; })) {
        r_print("no errors detected");
        return;
    }

    // init evaluator_poly with synds
    auto evaluator_poly = &synds[0];  // size = ecc_len * 2
    auto locator_poly = &temp_ntt1_ecc6[ecc_len * 2];  // size = ecc_len * 2
    _error_locator(&error_pos_rbo[0], error_count, &locator_poly[0]);
    r_print_vec("locator_poly", locator_poly, ecc_len);

    auto locator_poly_deg = error_count;
    r_print_vec("locator_poly_deg", &locator_poly_deg, 1);

    _error_evaluator(&locator_poly[0], locator_poly_deg, &evaluator_poly[0], &temp_ntt1_ecc6[ecc_len * 4]);
    r_print_vec("evaluator_poly", evaluator_poly, ecc_len);

    auto evaluator_poly_len = ecc_len;
    auto locator_poly_deriv = &temp_ntt1_ecc6[ecc_len * 1];

    // TODO: remove copy_n and add dst array to _deriv
    std::copy_n(&locator_poly[0], ecc_len, &locator_poly_deriv[0]);
    _deriv_shifted(locator_poly_deriv, ecc_len);
    auto locator_poly_deriv_len = _get_vec_size(&locator_poly_deriv[0], ecc_len);
    r_print_vec("locator_poly_deriv", locator_poly_deriv, ecc_len);

    // TODO: only need to check one lane
    auto locator_poly_deriv_len_max = vec<GFT>::max(locator_poly_deriv_len);

    for (size_t i = 0; i < error_count; ++i) {
        auto pos_rbo = error_pos_rbo[i];
        auto pos = ntt.rbo(pos_rbo);
        r_print("error pos:", pos);

        ::GFT x = ntt._roots_ntt[pos_rbo];

        GFT error = _forney(
            &locator_poly_deriv[0], locator_poly_deriv_len_max,
            &evaluator_poly[0], evaluator_poly_len,
            x
        );

        r_print_vec("error", &error, 1);

        block[pos] = gf.add(block[pos], error);
    }
    r_print("repair <");
}


template<size_t W>
void RSi16v<W>::_repair_ntt(GFT *const block, const size_t *const error_pos_rbo, size_t error_count, GFT *const temp_ntt1_ecc6) const {
        ntt.ntt(block);
        auto locator_poly = &temp_ntt1_ecc6[block_len];
        _error_locator_rev(error_pos_rbo, error_count, locator_poly);

        auto err_ntt = &temp_ntt1_ecc6[0];

        // copy synds
        std::copy_n(&block[0], ecc_len, &err_ntt[0]);

        for (size_t j = ecc_len; j < ntt_len; ++j) {
            GFT sum = GFT{0};
            for (size_t i = 0; i < error_count; ++i)
                sum = gf.sub(sum, gf.mul(locator_poly[i], err_ntt[j - error_count + i]));

            err_ntt[j] = sum;
            block[j] = gf.sub(block[j], sum);
        }

        std::fill_n(&block[0], ecc_len, GFT{0});

        // TODO limit to error positions only, test if performance improvement
        // for (size_t j = 0; j < error_count; ++j) {
        //     size_t i = error_pos_rbo[j];
        //     block[i] = gf.div(block[i], ntt_len);
        // }
        ntt.intt(block);
}

template<size_t W>
RSi16v<W>::GFT RSi16v<W>::_sugiyama(GFT *const a1, GFT *const r1, GFT *const temp_ecc4) const {
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
    // ntt.nttr(a2);
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
        auto r1h = vec<GFT>::gather(&r1[0], (r1_len - 1) & ::GFT(ecc_len_mask));
        auto r1h2 = vec<GFT>::gather(&r1[0], GFT(r1_len - 2) & ::GFT(ecc_len_mask));

        a1[1] = gf.neg(gf.inv(r1h));
        a1[0] = gf.mul(r1h2, gf.mul(a1[1], a1[1]));
        a1_len = GFT{} + 2;

        for (size_t i = 0; i < ecc_len; ++i)
            r1[i] = gf.mul(a1[0], r1[i]);
        for (size_t i = 1; i < ecc_len; ++i)
            r1[i] = gf.add(r1[i], gf.mul(a1[1], r2[i - 1]));

        r1_len = _get_vec_size(r1, ecc_len - 1);
        vec<GFT>::fill_n(&r1[0], r1_len, ::GFT(ecc_len) - r1_len, GFT{0});

        ld_print_vec("r1_len", &r1_len, 1);
        ld_print_vec("r1", r1, ecc_len);

        ld_print_vec("a1", a1, 2);
        ntt.nttr(a1);
        ntt.nttr(r2);
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
    if (vec<GFT>::max(r1_len) >= ecc_len / 2) {
        for (size_t i = 1; i < ecc_len / 2; ++i) {
            ld_print("------\nIteration begin");

            // auto r1h = r1[r1_len - 1];
            // TODO: mask gather
            auto r1h = vec<GFT>::gather(&r1[0], (r1_len - 1) & ::GFT(ecc_len_mask));
            auto r1l = r1[0];
            ntt.nttr(r1);

            ld_print_vec_intt("r2", r2, ecc_len);
            ld_print_vec_intt("r1", r1, ecc_len);
            ld_print_vec_intt("a2", a2, ecc_len);
            ld_print_vec_intt("a1", a1, ecc_len);

            // q = r2 // r1
            q_len = ntt.poly_div_ntt(r2, r2_len, r1, r1_len, r2l, r1h, t, q);

            ld_print("------");
            ld_print_vec("q_len", &q_len, 1);
            ld_print_vec_intt("q", q, ecc_len);

            // t = q * a1
            t_len = ntt.poly_mul_ntt(q, q_len, a1, a1_len, t);
            // q = q * r1
            q_len = ntt.poly_mul_ntt(q, q_len, r1, r1_len, q);

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

                if (vec<GFT>::is_zero(cond)) {
                    ntt.inttr(r1);
                    break;
                }

                // a1 = t;
                // r1 = q;
                vec<GFT>::copy_n_masked(&t[0], ecc_len, &a1[0], cond);
                vec<GFT>::copy_n_masked(&q[0], ecc_len, &r1[0], cond);
                a1_len = t_len;
                r1_len = q_len;

                if constexpr (!std::is_integral_v<GFT>)
                    lanes &= cond;
            }

            ntt.inttr(r1);

            r1_len = _get_vec_size(r1, vec<GFT>::max(r1_len & lanes));
            ld_print_vec("r1_len", &r1_len, 1);
            ld_print("------\nIteration end");
        }
    }

    ntt.inttr(a1);

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


template<size_t W>
inline void RSi16v<W>::_find_roots_ntt(
    GFT *const block,
    const GFT *const locator_poly, GFT const& locator_poly_len,
    const GFT *const locator_poly_deriv, GFT const& locator_poly_deriv_len,
    const GFT *const evaluator_poly, GFT const& evaluator_poly_len,
    GFT *const roots
) const {
    std::fill_n(&roots[0], ntt_len, GFT{0});
    std::copy_n(&locator_poly[0], ecc_len, &roots[0]);

    py_assert(vec<GFT>::max(locator_poly_len) <= ecc_len);
    ntt.pintt_ecc(&roots[0]);

    GFT lanes = ~GFT{0};
    if constexpr (!std::is_integral_v<GFT>)
        lanes = (locator_poly_len != 0);

    auto locator_poly_deriv_len_max = vec<GFT>::max(locator_poly_deriv_len & lanes);
    auto evaluator_poly_len_max = vec<GFT>::max(evaluator_poly_len & lanes);

    r_print("locator_poly_deriv_len_max", locator_poly_deriv_len_max);
    r_print("evaluator_poly_len_max", evaluator_poly_len_max);

    for (::GFT i = 0; i < block_len; ++i) {
        auto i_rbo = ntt.rbo(i);
        // auto x = gf.pow(root, i_rbo);
        ::GFT x = ntt._roots_ntt[i_rbo];
        GFT root_locations = (roots[i] == 0) & lanes;

        if (vec<GFT>::any(root_locations)) {
            GFT error = _forney(
                &locator_poly_deriv[0], locator_poly_deriv_len_max,
                &evaluator_poly[0], evaluator_poly_len_max,
                x
            );

            if constexpr (!std::is_integral_v<GFT>)
                error &= root_locations;

            r_print("i:", i, "rbo(i):", i_rbo);
            r_print_vec("error", &error, 1);

            block[i] = gf.add(block[i], error);
        }
    }
}


template<size_t W>
RSi16v<W>::GFT RSi16v<W>::_forney(
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


template<size_t W>
inline void RSi16v<W>::_error_locator_rev(const size_t *const error_pos_rbo, size_t error_count, GFT *const locator_poly) const {
    std::fill_n(&locator_poly[0], ecc_len, GFT{0});

    size_t last_error = error_count - 1;

    for (size_t j = 0; j < error_count; ++j) {
        GFT x = GFT{} + gf.pow(root, error_pos_rbo[j]);
        for (size_t i = last_error - j; i < last_error; ++i) {
            locator_poly[i] = gf.sub(locator_poly[i], gf.mul(locator_poly[i + 1], x));
        }
        locator_poly[last_error] = gf.sub(locator_poly[last_error], x);
    }
}


template<size_t W>
inline void RSi16v<W>::_error_locator(const size_t *const error_pos_rbo, size_t error_count, GFT *const locator_poly) const {
    std::fill_n(&locator_poly[0], ecc_len, GFT{0});

    for (size_t j = 0; j < error_count; ++j) {
        GFT x = GFT{} + ntt._roots_ntt[error_pos_rbo[j]];
        for (size_t i = error_count - 1; i > 0; --i) {
            locator_poly[i] = gf.sub(locator_poly[i], gf.mul(locator_poly[i - 1], x));
        }
        locator_poly[0] = gf.sub(locator_poly[0], x);
    }
}


template<size_t W>
inline void RSi16v<W>::_error_evaluator(const GFT *const locator_poly, size_t locator_poly_deg, GFT *const evaluator_poly, GFT *const temp) const {
    // evaluator_poly initialized with synds
    std::fill_n(&evaluator_poly[ecc_len], ecc_len, GFT{0});

    ntt.nttr2(&evaluator_poly[0]);

    std::copy_n(&locator_poly[0], locator_poly_deg, &temp[1]);
    temp[0] = GFT{} + 1;
    std::fill_n(&temp[locator_poly_deg + 1], ecc_len * 2 - locator_poly_deg, GFT{0});

    ntt.nttr2(&temp[0]);

    for (size_t i = 0; i < ecc_len * 2; ++i)
        evaluator_poly[i] = gf.mul(evaluator_poly[i], temp[i]);

    ntt.inttr2(&evaluator_poly[0]);
    std::fill_n(&evaluator_poly[ecc_len], ecc_len, GFT{0});
}


template<size_t W>
inline size_t RSi16v<W>::_vec_shift(const GFT *const a, size_t a_len, size_t shift, GFT *const r) const {
    // ::GFT ecc_root = gf.pow(root, ntt_len / ecc_len);

    for (size_t i = 0; i < ecc_len; ++i)
        // r[i] = gf.mul(a[i], gf.pow(ecc_root, shift * ntt._rbo_ecc[i] & ecc_len_mask));
        r[i] = gf.mul(a[i], ntt._roots_ecc[shift * ntt._rbo_ecc[i] & ecc_len_mask]);

    return a_len + shift;
}


template<size_t W>
inline RSi16v<W>::GFT RSi16v<W>::_vec_sub(const GFT *const a, GFT a_len, GFT *const b, GFT b_len, GFT *const r) const {
    for (size_t i = 0; i < ecc_len; ++i)
        r[i] = gf.sub(a[i], b[i]);

    return vec<GFT>::max(a_len, b_len);
}


template<size_t W>
inline size_t RSi16v<W>::_norm_size_rev(GFT *const r, size_t r_len) const {
    check_le_ecc(r_len);
    size_t start = r_len;
    while (start < ecc_len && vec<GFT>::is_zero(r[start]))
        ++start;

    if (start == ecc_len) {
        start = 0;
        while (start < r_len && vec<GFT>::is_zero(r[start]))
            ++start;
    }

    std::rotate(&r[0], &r[start], &r[ecc_len]);
    if (r_len == ecc_len && start == 0)
        return r_len;
    else
        return (r_len - start) & ecc_len_mask;
}


template<size_t W>
inline RSi16v<W>::GFT RSi16v<W>::_get_vec_size(const GFT *const r, size_t r_len_max) const {
    check_le_ecc(r_len_max);

    GFT len = GFT{} + ::GFT(r_len_max);

    if constexpr (std::is_integral_v<GFT>) {
        while (r_len_max > 0 && vec<GFT>::is_zero(r[r_len_max - 1]))
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


template<size_t W>
inline size_t RSi16v<W>::_deriv(GFT *const r, size_t r_len) const {
    check_le_ecc(r_len);

    for (::GFT i = 0; i < r_len - 1; ++i)
        r[i] = gf.mul(r[i + 1], i + 1);
    r[r_len - 1] = GFT{0};
    return r_len - 1;
}


template<size_t W>
inline size_t RSi16v<W>::_deriv_shifted(GFT *const r, size_t r_len) const {
    check_le_ecc(r_len);

    for (::GFT i = 0; i < r_len; ++i)
        r[i] = gf.mul(r[i], i + 1);
    return r_len;
}

template<size_t W>template<typename T>
inline T RSi16v<W>::_eval(const T *const r, size_t r_len, T x) const {
    check_le_ecc(r_len);

    T sum = T{0};
    for (size_t i = r_len - 1; i < r_len; --i)
        sum = gf.add(gf.mul(sum, x), r[i]);
    return sum;
}


template<size_t W>
inline void RSi16v<W>::_reverse_ntt(GFT *const vec, size_t shift) const {
    check_le_ecc(shift);

    for (size_t i = 1; i < ecc_len; ++i) {
        // ntt.rbo(ecc_len - ntt.rbo(i))
        auto j = ntt._rbo_ecc[i];
        auto k = ntt._rbo_ecc[ecc_len - j];
        if (i < k)
            std::swap(vec[i], vec[k]);
        vec[i] = gf.mul(vec[i], ntt._roots_i_ecc[j]);
        vec[i] = gf.mul(vec[i], ntt._roots_ecc[shift * j & ecc_len_mask]);
    }
}


template<size_t W>
inline void RSi16v<W>::_reverse_ntt_vec(GFT *const vec, GFT shift) const {
    for (size_t i = 1; i < ecc_len; ++i) {
        // ntt.rbo(ecc_len - ntt.rbo(i))
        auto j = ntt._rbo_ecc[i];
        auto k = ntt._rbo_ecc[ecc_len - j];
        if (i < k)
            std::swap(vec[i], vec[k]);
        vec[i] = gf.mul(vec[i], ntt._roots_i_ecc[j]);
        vec[i] = gf.mul(vec[i], gf.gather(&ntt._roots_ecc[0], shift * j & ::GFT(ecc_len_mask)));
    }
}
