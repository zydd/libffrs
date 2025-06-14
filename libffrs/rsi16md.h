/**************************************************************************
 * rsi16md.h
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

# pragma once

#include <array>
#include <vector>

template<size_t N>
class RSi16v {
public:
    RSi16v(size_t block_size, size_t ecc_len, uint32_t primitive);
    void encode(uint32_t block[]) const;

protected:
    struct data;
    data *d;
};


// class RSi16v16 {
// public:
//     RSi16v16(size_t block_size, size_t ecc_len, uint32_t primitive);
//     void encode(uint32_t block[]) const;

// protected:
//     struct data;
//     data *d;
// };
