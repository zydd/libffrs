/**************************************************************************
 * rsi16v.hpp
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

#include <array>
#include <cstdint>
#include <vector>

#include "pygfi16.hpp"

using GFT = uint32_t;

#include "ntt.hpp"


template<size_t W>
class RSi16v {
private:
    const GFi16 &gf;
    const NTT ntt;
    std::vector<::GFT> _ecc_mix;
    std::vector<::GFT> _ecc_mix_i;

public:
    const ::GFT root;
    const size_t ntt_len;
    const size_t block_len;
    const size_t ecc_len;
    const ::GFT ecc_len_mask;
    using GFT = simd_map_t<W>;
    using cGFT = const GFT *const;
    using mGFT = GFT *const;

public:
    RSi16v(GFi16 const& gf, size_t block_len, size_t ecc_len);
    ~RSi16v();

    inline void encode(::GFT block[]) const
        { _encode(reinterpret_cast<GFT *>(block)); }

    inline void repair(::GFT block[], const size_t error_pos_rbo[], size_t error_count, ::GFT temp_ntt1_ecc6[]) const
        { _repair(reinterpret_cast<GFT *>(block), error_pos_rbo, error_count, reinterpret_cast<GFT *>(temp_ntt1_ecc6)); }

    inline void repair_ntt(::GFT block[], const size_t error_pos_rbo[], size_t error_count, ::GFT temp_ntt1_ecc6[]) const
        { _repair_ntt(reinterpret_cast<GFT *>(block), error_pos_rbo, error_count, reinterpret_cast<GFT *>(temp_ntt1_ecc6)); }

    inline void repair(::GFT block[], ::GFT temp_ntt1_ecc6[]) const
        { _repair(reinterpret_cast<GFT *>(block), reinterpret_cast<GFT *>(temp_ntt1_ecc6)); }

    inline void pntt(::GFT block[]) const
        { ntt.pntt(reinterpret_cast<GFT *>(block)); }

    inline void mix_ecc(::GFT ecc[]) const
        { _mix_ecc(reinterpret_cast<GFT *>(ecc)); }

    uint16_t rbo(uint16_t) const;

protected:
    void _encode(GFT block[]) const;
    void _repair(GFT block[], GFT temp_ntt1_ecc6[]) const;
    void _repair(GFT block[], const size_t error_pos_rbo[], size_t error_count, GFT temp_ntt1_ecc6[]) const;
    void _repair_ntt(GFT block[], const size_t error_pos_rbo[], size_t error_count, GFT temp_ntt1_ecc6[]) const;
    void _mix_ecc(GFT ecc[]) const;
    void _mix_ecc_residue(GFT ecc[]) const;

private:
    GFT _sugiyama(GFT *a1, GFT *r1, GFT *temp_ecc4) const;
    GFT _forney(const GFT *locator_poly_deriv, size_t locator_poly_deriv_len, const GFT *evaluator_poly, size_t evaluator_poly_len, const ::GFT& x) const;
    void _error_locator_rev(const size_t *error_pos_rbo, size_t error_count, GFT *locator_poly) const;
    void _error_locator(const size_t *error_pos_rbo, size_t error_count, GFT *locator_poly) const;
    void _error_evaluator(const GFT *locator_poly, size_t locator_poly_deg, GFT *evaluator_poly, GFT *temp) const;
    void _find_roots_ntt(GFT *block, const GFT *locator_poly, const GFT& locator_poly_len, const GFT *locator_poly_deriv, const GFT& locator_poly_deriv_len, const GFT *evaluator_poly, const GFT& evaluator_poly_len, GFT *roots) const;
    size_t _vec_shift(const GFT *a, size_t a_len, size_t shift, GFT *r) const;
    GFT _vec_sub(const GFT *a, GFT a_len, GFT *b, GFT b_len, GFT *r) const;
    size_t _norm_size_rev(GFT *r, size_t r_len) const;
    GFT _get_vec_size(const GFT *r, size_t r_len_max) const;
    size_t _deriv(GFT *r, size_t r_len) const;
    size_t _deriv_shifted(GFT *r, size_t r_len) const;

    template<typename T>
    T _eval(const T* r, size_t r_len, T x) const;

    void _reverse_ntt(GFT *vec, size_t shift) const;
    void _reverse_ntt_vec(GFT *vec, GFT shift) const;
};
