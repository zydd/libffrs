/**************************************************************************
 * rsi16md.cpp
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


#include "rsi16md_impl.hpp"


template class RSi16v<1>;


template<>
simd_mask_t RSi16vImpl<GFT>::is_zero(GFT const& vec) {
    return (vec == 0);
}


template<>
void RSi16vImpl<GFT>::assign_masked(GFT& vec, GFT const& value, GFT const& condition) const {
    if (condition)
        vec = value;
}


template<>
void RSi16vImpl<GFT>::copy_n_masked(const GFT src[], size_t n, GFT dst[], GFT const& condition) const {
    for (size_t i = 0; i < n; i++) {
        if (condition)
            dst[i] = src[i];
    }
}
