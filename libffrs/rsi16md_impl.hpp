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

    inline RSi16vImpl(GFi16 const& gf, size_t block_len, size_t ecc_len):
        gf(gf),
        block_len(block_len),
        ecc_len(ecc_len)
    {
        _rbo_shift = 16 - __builtin_ctzl(block_len);
        block_len_iv = GFT{0} + gf.inv(block_len);

        root = gf.exp(gf.div(gf.log(1), block_len));
        if (root >= 0x8000)
            root = gf.neg(root);

        if (gf.pow(root, block_len) != 1)
            throw std::runtime_error("Root of unity not found for block size");

        uint32_t root_i = gf.inv(root);
        _roots_v_block.resize(block_len);
        _roots_iv_block.resize(block_len);
        for (size_t i = 0; i < block_len; ++i) {
            _roots_v_block[i] = GFT{0} + gf.pow(root, i);
            _roots_iv_block[i] = GFT{0} + gf.pow(root_i, i);
        }

        _ecc_mix_v.resize(ecc_len);
        for (size_t i = 0; i < ecc_len; ++i) {
            auto w_i = gf.inv(*reinterpret_cast<uint32_t *>(&_roots_v_block[rbo(block_len - ecc_len)]));
            _ecc_mix_v[i] = GFT{0} + gf.neg(gf.div(gf.pow(w_i, i), ecc_len));
        }

        // uint32_t ecc_root = gf.inv(gf.pow(root, block_len / ecc_len));
        uint32_t ecc_root = gf.inv(gf.exp(gf.div(gf.log(1), ecc_len)));
        _roots_iv_ecc.resize(ecc_len);
        for (size_t i = 0; i < ecc_len; ++i) {
            _roots_iv_ecc[i] = GFT{0} + gf.pow(ecc_root, i);
        }
    }

    inline void encode(GFT block[]) const {
        ct_butterfly(&_roots_v_block[0], &block[0], block_len);

        for (size_t j = 0; j < ecc_len; ++j)
            block[j] = gf.mul(block[j], _ecc_mix_v[j]);

        gs_butterfly(&_roots_iv_ecc[0], &block[0], ecc_len, ecc_len);
    }

    inline void ntt(GFT block[]) const {
        ct_butterfly(&_roots_v_block[0], &block[0], block_len);
    }

    inline void repair(GFT block[], const size_t error_pos_rbo[], size_t error_count, GFT temp[]) const {
        ct_butterfly(&_roots_v_block[0], &block[0], block_len);

        // copy synds
        std::copy_n(&block[0], ecc_len, &temp[0]);

        // re-use start of block to store error locator polynomial (omiting 1)
        // deg(locator_poly) == error_count <= ecc_len
        auto locator_poly = block;
        std::fill_n(&locator_poly[0], ecc_len, GFT{0});

        size_t last_error = error_count - 1;

        for (size_t j = 0; j < error_count; ++j) {
            GFT x = GFT{0} + gf.neg(gf.pow(root, error_pos_rbo[j]));
            for (size_t i = last_error - j; i < last_error; ++i) {
                locator_poly[i] = gf.add(locator_poly[i], gf.mul(locator_poly[i + 1], x));
            }
            locator_poly[last_error] = gf.add(locator_poly[last_error], x);
        }

        for (size_t j = ecc_len; j < block_len; ++j) {
            GFT sum = GFT{0};
            for (size_t i = 0; i < error_count; ++i)
                sum = gf.sub(sum, gf.mul(locator_poly[i], temp[j - error_count + i]));

            temp[j] = sum;
            block[j] = gf.sub(block[j], sum);
        }

        std::fill_n(&block[0], ecc_len, GFT{0});
        gs_butterfly(&_roots_iv_block[0], &block[0], block_len, block_len);
        for (size_t j = 0; j < block_len; ++j)
            block[j] = gf.mul(block[j], block_len_iv);

        // TODO limit to error positions only
        // for (size_t j = 0; j < error_count; ++j) {
        //     size_t i = error_pos_rbo[j];
        //     block[i] = gf.div(block[i], block_len);
        // }
    }

    uint16_t rbo(uint16_t b) const {
        return ffrs::detail::rbo16(b) >> _rbo_shift;
    }

protected:
    GFi16 const& gf;
    uint16_t _rbo_shift;
    size_t block_len;
    size_t ecc_len;
    GFT block_len_iv;
    std::vector<GFT> _roots_v_block;
    std::vector<GFT> _roots_iv_block;
    std::vector<GFT> _roots_iv_ecc;
    std::vector<GFT> _ecc_mix_v;

    inline void ct_butterfly(const GFT roots[], GFT block[], size_t block_len) const {
        for (size_t stride = 1, exp_f = block_len >> 1; stride < block_len; stride *= 2, exp_f >>= 1) {
            for (size_t start = 0; start < block_len /*input_size*/; start += stride * 2) {
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

    inline void gs_butterfly(const GFT roots[], GFT block[], size_t block_len, size_t end) const {
        for (size_t stride = block_len / 2, exp_f = 0; stride > 0; stride /= 2, exp_f += 1) {
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
};
