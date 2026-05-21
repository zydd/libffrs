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


template<>
void RSi16vImpl<GFTx8>::add_masked(GFTx8 vec[], GFTx8 const& i, GFTx8 const& value, GFTx8 const& condition) const {
    auto index = i * 8 + GFTx8{0, 1, 2, 3, 4, 5, 6, 7};
    auto prev_value = (GFTx8) _mm256_mask_i32gather_epi32(
        __m256i{},
        (const int *) vec,
        (__m256i) index,
        (__m256i) condition,
        4
    );

    auto sum = gf.add(prev_value, value);

    for (int j = 0; j < 8; j++)
        if (condition[j])
            vec[i[j]][j] = sum[j];
}


template<>
simd_mask_t RSi16vImpl<GFTx8>::is_zero(GFTx8 const& vec) {
    // return std::all_of(
    //     reinterpret_cast<const uint32_t *>(&vec),
    //     reinterpret_cast<const uint32_t *>(&vec) + 8,
    //     [](auto v) { return v == 0; }
    // );
    return _mm256_testz_si256((__m256i) vec, (__m256i) vec);
}
