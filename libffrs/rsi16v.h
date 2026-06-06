/**************************************************************************
 * rsi16v.h
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

# pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "pygfi16.hpp"

using GFT = uint32_t;


template<size_t W>
class RSi16v {
public:
    RSi16v(GFi16 const& gf, size_t block_size, size_t ecc_len);
    ~RSi16v();
    void encode(GFT block[]) const;
    void repair(GFT block[], const size_t error_pos_rbo[], size_t error_count, GFT temp2[]) const;
    void repair_ntt(GFT block[], const size_t error_pos_rbo[], size_t error_count, GFT temp2[]) const;
    void repair(GFT block[], GFT temp2[]) const;
    void pntt(GFT block[]) const;
    void ecc_mix(GFT block[]) const;
    uint16_t rbo(uint16_t) const;

protected:
    class Impl;
    Impl *d;

public:
    const GFT root;
    const size_t ntt_len;
};
