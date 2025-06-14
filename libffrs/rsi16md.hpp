/**************************************************************************
 * rsi16md.hpp
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

template<typename GF>
class RSi16v {
public:
    using GFT = typename GF::GFT;

    RSi16v(size_t block_size, size_t ecc_len, uint32_t primitive);

    void encode(GFT block[], size_t block_size) const;

protected:
    GF gf;
    size_t block_size;
    size_t ecc_len;
    std::vector<GFT> _rootsv;
    std::vector<GFT> _roots_iv;
    std::vector<uint16_t> _rbo;
    std::vector<GFT> _ecc_mix_wv;

    void ct_butterfly(const GFT roots[], GFT block[], size_t block_size) const;

    void gs_butterfly(size_t end, const GFT roots[], GFT block[], size_t block_size) const;
};
