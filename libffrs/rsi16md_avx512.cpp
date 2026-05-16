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
GFTx16 RSi16vImpl<GFTx16>::gather(const uint32_t vec[], GFTx16 const& i) const {
    return (GFTx16) _mm512_i32gather_epi32(
        (__m512i) i,
        vec,
        4
    );
}
