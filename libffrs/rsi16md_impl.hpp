/**************************************************************************
 * rsi16md_impl.hpp
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

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

#include "pygfi16.hpp"


template<typename GFT>
class RSi16vImpl {
public:
    uint32_t root;

    inline RSi16vImpl(GFi16 const& gf, size_t block_size, size_t ecc_len):
        gf(gf),
        block_size(block_size),
        ecc_len(ecc_len)
    {
        size_t nbits = __builtin_ctzl(block_size);

        root = gf.exp(gf.div(gf.log(1), block_size));

        if (root >= 0x8000)
            root = gf.neg(root);

        if (gf.pow(root, block_size) != 1)
            throw std::runtime_error("Root of unity not found for block size");

        _rootsv.resize(block_size);
        _roots_iv.resize(block_size);
        uint32_t r = root;
        uint32_t r_i = gf.inv(r);
        for (size_t i = 0; i < block_size; ++i) {
            // _roots[i] = r;
            // r = gf.pow(r, 2);
            _rootsv[i] = GFT{0} + gf.pow(root, i);
            _roots_iv[i] = GFT{0} + gf.pow(r_i, i);
        }

        _rbo.resize(block_size);
        for (size_t i = 0; i < block_size; ++i)
            _rbo[i] = rbo16(i) >> (16 - nbits);

        _ecc_mix_wv.resize(ecc_len);
        for (size_t j = 0; j < ecc_len; ++j) {
            auto w = *reinterpret_cast<uint32_t *>(&_roots_iv[_rbo[block_size - ecc_len]]);
            // w = - w ** (- rbo(i) * j) / block_size
            _ecc_mix_wv[j] = GFT{0} + gf.neg(gf.div(gf.pow(w, j), block_size));
        }
    }

    inline void encode(GFT block[]) const {
        ct_butterfly(&_rootsv[0], block, block_size);

        for (size_t j = 0; j < ecc_len; ++j)
            block[j] = gf.mul(block[j], _ecc_mix_wv[j]);

        for (size_t j = 1; j < block_size / ecc_len; ++j)
            std::copy_n(&block[0], ecc_len, &block[j * ecc_len]);
            // std::memcpy(&block[j * ecc_len], block, ecc_len * sizeof(GFT));

        gs_butterfly(ecc_len, &_roots_iv[0], block, block_size);
    }

protected:
    GFi16 const& gf;
    size_t block_size;
    size_t ecc_len;
    std::vector<GFT> _rootsv;
    std::vector<GFT> _roots_iv;
    std::vector<uint16_t> _rbo;
    std::vector<GFT> _ecc_mix_wv;

    inline void ct_butterfly(const GFT roots[], GFT block[], size_t block_size) const {
        for (size_t stride = 1, exp_f = block_size >> 1; stride < block_size; stride *= 2, exp_f >>= 1) {
            for (size_t start = 0; start < block_size /*input_size*/; start += stride * 2) {
                {
                    // Cooley-Tukey butterfly
                    GFT a = block[start];
                    GFT b = block[start + stride];
                    block[start] = gf.add(a, b);
                    block[start + stride] = gf.sub(a, b);
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // j = i - start
                    GFT w = roots[exp_f * (i - start)];

                    // Cooley-Tukey butterfly
                    GFT a = block[i];
                    GFT b = block[i + stride];
                    GFT m = gf.mul(b, w);
                    block[i] = gf.add(a, m);
                    block[i + stride] = gf.sub(a, m);
                }
            }
        }
    }

    inline void gs_butterfly(size_t end, const GFT roots[], GFT block[], size_t block_size) const {
        for (size_t stride = block_size / 2, exp_f = 0; stride > 0; stride /= 2, exp_f += 1) {
            for (size_t start = 0; start < end; start += stride * 2) {
                {
                    // Gentleman-Sande butterfly
                    GFT a = block[start];
                    GFT b = block[start + stride];
                    block[start] = gf.add(a, b);
                    block[start + stride] = gf.sub(a, b);
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // Gentleman-Sande butterfly
                    GFT w = roots[(i - start) << exp_f];

                    GFT a = block[i];
                    GFT b = block[i + stride];
                    block[i] = gf.add(a, b);
                    block[i + stride] = gf.mul(gf.sub(a, b), w);
                }
            }
        }
    }

    uint16_t rbo16(uint16_t b) const {
        b = ((b >> 1) & 0x5555) | ((b & 0x5555) << 1);
        b = ((b >> 2) & 0x3333) | ((b & 0x3333) << 2);
        b = ((b >> 4) & 0x0f0f) | ((b & 0x0f0f) << 4);
        b = ((b >> 8)         ) | ((b         ) << 8);
        return b;
    }
};
