/**************************************************************************
 * rsi16md_avx512.cpp
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

template class RSi16v<16>;


template<>
GFTx16 ffrs::simd_gather_base::gather(const uint32_t vec[], GFTx16 const& i) {
    return (GFTx16) _mm512_i32gather_epi32(
        (__m512i) i,
        vec,
        4
    );
}


template<>
void RSi16vImpl<GFTx16>::add_masked(GFTx16 vec[], GFTx16 const& i, GFTx16 const& value, GFTx16 const& condition) const {
    auto mask = _mm512_movepi32_mask((__m512i) condition);
    auto index = i * 16 + GFTx16{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    auto prev_value = (GFTx16) _mm512_mask_i32gather_epi32(
        __m512i{},
        mask,
        (__m512i) index,
        vec,
        4
    );

    _mm512_mask_i32scatter_epi32(
        vec,
        mask,
        (__m512i) index,
        (__m512i) gf.add(prev_value, value),
        4
    );
}
