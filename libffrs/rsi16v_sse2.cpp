/**************************************************************************
 * rsi16v_sse2.cpp
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

#include <cstdint>
#include <immintrin.h>

#include "rsi16v_impl.hpp"
#include "simd.hpp"


template class RSi16v<4>;


template<>
GFTx4 ffrs::simd_gather_base::gather(const uint32_t vec[], GFTx4 const& i) {
    // Requires AVX-2
    // return (GFTx4) _mm_i32gather_epi32(
    //     (const int*) vec,
    //     (__m128i) i,
    //     4
    // );
    return GFTx4{
        vec[i[0]],
        vec[i[1]],
        vec[i[2]],
        vec[i[3]]
    };
}


template<>
void vec<GFTx4>::assign_masked(GFTx4& vec, GFTx4 const& value, GFTx4 const& condition) {
    for (int j = 0; j < 4; j++)
        if (condition[j])
            vec[j] = value[j];
}


template<>
void vec<GFTx4>::copy_n_masked(const GFTx4 src[], size_t n, GFTx4 dst[], GFTx4 const& condition) {
    for (size_t i = 0; i < n; i++) {
        for (int j = 0; j < 4; j++)
            if (condition[j])
                dst[i][j] = src[i][j];
    }
}


template<>
simd_mask_t vec<GFTx4>::is_zero(GFTx4 const& vec) {
    // return !(vec[0] || vec[1] || vec[2] || vec[3]);
    return _mm_movemask_epi8(_mm_cmpeq_epi32((__m128i) vec, _mm_setzero_si128())) == 0xffff;
}
