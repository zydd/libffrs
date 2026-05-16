/**************************************************************************
 * rsi16md_sse2.cpp
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

#include "rsi16md_impl.hpp"

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

// template<>
// GFTx4 ffrs::simd_gather_base::gather(const GFTx4 vec[], GFTx4 const& i) {
//     py_assert(!"GFTx4 ffrs::simd_gather_base::gather");
// }


// template<>
// void ffrs::simd_gather_base::scatter(GFTx4 vec[], GFTx4 const& i, GFTx4 const& value, GFTx4 const& mask) {
//     _mm_mask_i32scatter_epi32(
//         vec,
//         _mm_movemask_epi8((__m128i) mask),
//         (__m128i) i,
//         (__m128i) value,
//         4
//     );
// }

