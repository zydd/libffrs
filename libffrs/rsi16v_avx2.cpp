/**************************************************************************
 * rsi16v_avx2.cpp
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

#define RSI16V_IMPL_INSTANCE_W 8
#include "rsi16v_impl.hpp"

#include "simd.hpp"


template<>
GFTx8 ffrs::simd_gather_base::gather(const uint32_t vec[], GFTx8 const& i) {
    return (GFTx8) _mm256_i32gather_epi32(
        (const int *) vec,
        (__m256i) i,
        4
    );
}


// template<>
// void vec::assign_masked<GFTx8>(GFTx8& vec, GFTx8 const& value, GFTx8 const& condition) {
//     for (int j = 0; j < 8; j++)
//         if (condition[j])
//             vec[j] = value[j];
// }


template<>
void vec::copy_n_masked<GFTx8>(const GFTx8 src[], size_t n, GFTx8 dst[], GFTx8 const& condition) {
    // TODO: test if copying columns perform better
    for (size_t i = 0; i < n; i++) {
        for (int j = 0; j < 8; j++)
            if (condition[j])
                dst[i][j] = src[i][j];
    }
}


template<>
simd_mask_t vec::is_zero<GFTx8>(GFTx8 const& vec) {
    // return std::all_of(
    //     reinterpret_cast<const uint32_t *>(&vec),
    //     reinterpret_cast<const uint32_t *>(&vec) + 8,
    //     [](auto v) { return v == 0; }
    // );
    return _mm256_testz_si256((__m256i) vec, (__m256i) vec);
}
