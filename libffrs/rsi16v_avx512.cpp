/**************************************************************************
 * rsi16v_avx512.cpp
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

#define RSI16V_IMPL_INSTANCE_W 16
#include "rsi16v_impl.hpp"

#include "simd.hpp"


template<>
GFTx16 ffrs::simd_gather_base::gather(const uint32_t vec[], GFTx16 const& i) {
    return (GFTx16) _mm512_i32gather_epi32(
        (__m512i) i,
        vec,
        4
    );
}


// template<>
// void vec::assign_masked<GFTx16>(GFTx16& vec, GFTx16 const& value, GFTx16 const& condition) {
//     for (int j = 0; j < 16; j++)
//         if (condition[j] || true)
//             vec[j] = value[j];

//     // auto mask = _mm512_movepi32_mask((__m512i) condition);
//     // vec = (GFTx16) _mm512_mask_mov_epi32(
//     //     (__m512i) vec,
//     //     mask,
//     //     (__m512i) value
//     // );
// }


template<>
void vec::copy_n_masked<GFTx16>(const GFTx16 src[], size_t n, GFTx16 dst[], GFTx16 const& condition) {
    // for (size_t i = 0; i < n; i++) {
    //     for (int j = 0; j < 16; j++)
    //         if (condition[j])
    //             dst[i][j] = src[i][j];
    // }
    auto mask = _mm512_movepi32_mask((__m512i) condition);
    for (size_t i = 0; i < n; i++) {
        dst[i] = (GFTx16) _mm512_mask_mov_epi32(
            (__m512i) dst[i],
            mask,
            (__m512i) src[i]
        );
    }
}


template<>
simd_mask_t vec::is_zero<GFTx16>(GFTx16 const& vec) {
    // return std::all_of(
    //     reinterpret_cast<const uint32_t *>(&vec),
    //     reinterpret_cast<const uint32_t *>(&vec) + 16,
    //     [](auto v) { return v == 0; }
    // );
    // return (_mm512_cmpeq_epi32_mask((__m512i) vec, __m512i{0}) == 0xffff);
    return _mm512_test_epi32_mask((__m512i) vec, (__m512i) vec) == 0;
}
