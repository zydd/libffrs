/**************************************************************************
 * rsi16md_avx512.cpp
 *
 * Copyright 2025 Gabriel Machado
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
#include "rsi16md.h"

typedef uint32_t u32x16 __attribute__((vector_size(16 * sizeof(uint32_t))));

template <size_t N>
struct RSi16v<N>::data : RSi16vImpl<GFi16v<u32x16>> { using RSi16vImpl::RSi16vImpl; };


template <size_t N>
RSi16v<N>::RSi16v(size_t block_size, size_t ecc_len, uint32_t primitive) {
    d = new data(block_size, ecc_len, primitive);
}


template <size_t N>
void RSi16v<N>::encode(uint32_t block[]) const {
    d->encode(reinterpret_cast<typename data::GFT *>(block));
}


template class RSi16v<16>;
