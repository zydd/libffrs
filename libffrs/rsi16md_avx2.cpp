/**************************************************************************
 * rsi16md_avx2.cpp
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


#include "galois.hpp"


#include "rsi16md_impl.hpp"

template class RSi16v<8>;


template<>
GFTx8 ffrs::simd_gather_base::gather(const uint32_t vec[], GFTx8 const& i) {
    return (GFTx8) _mm256_i32gather_epi32(
        (const int *) vec,
        (__m256i) i,
        4
    );
}

// template<>
// GFTx8 ffrs::simd_gather_base::gather(const GFTx8 vec[], GFTx8 const& i) {
//     py_assert(!"GFTx8 ffrs::simd_gather_base::gather");
// }


// template<>
// void ffrs::simd_gather_base::scatter(GFTx8 vec[], GFTx8 const& i, GFTx8 const& value, GFTx8 const& mask) {
//     _mm256_mask_i32scatter_epi32(
//         vec,
//         _mm256_movemask_epi8((__m256i) mask),
//         (__m256i) i,
//         (__m256i) value,
//         4
//     );
// }
