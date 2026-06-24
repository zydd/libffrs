/**************************************************************************
 * rsi16v.cpp
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


#define RSI16V_IMPL_INSTANCE_W 1
#include "rsi16v_impl.hpp"

#include "simd.hpp"


template<>
simd_mask_t vec::is_zero<GFT>(GFT const& vec) {
    return (vec == 0);
}


template<>
void vec::assign_masked<GFT>(GFT& vec, GFT const& value, GFT const& condition) {
    if (condition)
        vec = value;
}


template<>
void vec::copy_n_masked(const GFT src[], size_t n, GFT dst[], GFT const& condition) {
    for (size_t i = 0; i < n; i++) {
        if (condition)
            dst[i] = src[i];
    }
}
