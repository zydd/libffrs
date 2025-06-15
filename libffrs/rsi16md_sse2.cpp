/**************************************************************************
 * rsi16md_sse2.cpp
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


typedef uint32_t u32x4 __attribute__((vector_size(4 * sizeof(uint32_t))));

template <size_t N>
struct RSi16v<N>::data : RSi16vImpl<GFi16v<u32x4>> { using RSi16vImpl::RSi16vImpl; };


template <size_t N>
RSi16v<N>::RSi16v(GFi16 const& gf, size_t block_size, size_t ecc_len) {
    d = new data(gf, block_size, ecc_len);
}


template <size_t N>
void RSi16v<N>::encode(uint32_t block[]) const {
    d->encode(reinterpret_cast<GFi16v<u32x4> *>(block));
}


template class RSi16v<4>;
