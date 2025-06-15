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

template<typename T, size_t N>
struct vec {
    std::array<T, N> v;

#define _VEC_OP_BIN(op) \
    inline vec operator op(const vec& other) const { \
        vec res; \
        for (size_t i = 0; i < N; ++i) \
            res.v[i] = v[i] op other.v[i]; \
        return res; \
    } \
    template<typename V> \
    inline vec operator op(const V& other) const { \
        vec res; \
        for (size_t i = 0; i < N; ++i) \
            res.v[i] = v[i] op other; \
        return res; \
    } \
    template<typename V> \
    friend vec operator op(const V& lhs, const vec& rhs) { \
        vec res; \
        for (size_t i = 0; i < N; ++i) \
            res.v[i] = lhs op rhs.v[i]; \
        return res; \
    } \

    inline vec operator-() const {
        vec res;
        for (size_t i = 0; i < N; ++i)
            res.v[i] = -v[i];
        return res;
    }
    inline vec operator!() const {
        vec res;
        for (size_t i = 0; i < N; ++i)
            res.v[i] = !v[i];
        return res;
    }

    #define _VEC_OP_ASSIGN(op) \
        inline vec& operator op(const vec& other) { \
            for (size_t i = 0; i < N; ++i) \
                v[i] op other.v[i]; \
            return *this; \
        } \
        template<typename V> \
        inline vec& operator op(const V& other) { \
            for (size_t i = 0; i < N; ++i) \
                v[i] op other; \
            return *this; \
        } \

    _VEC_OP_BIN(+)
    _VEC_OP_BIN(-)
    _VEC_OP_BIN(*)
    _VEC_OP_BIN(/)
    _VEC_OP_BIN(>>)
    _VEC_OP_BIN(<<)
    _VEC_OP_BIN(&)
    _VEC_OP_BIN(|)
    _VEC_OP_BIN(^)
    _VEC_OP_BIN(&&)
    _VEC_OP_BIN(||)
    _VEC_OP_BIN(>)
    _VEC_OP_BIN(>=)
    _VEC_OP_BIN(<)
    _VEC_OP_BIN(<=)
    _VEC_OP_BIN(==)
    _VEC_OP_BIN(!=)

    _VEC_OP_ASSIGN(+=)
    _VEC_OP_ASSIGN(-=)
    _VEC_OP_ASSIGN(*=)
    _VEC_OP_ASSIGN(/=)
    _VEC_OP_ASSIGN(>>=)
    _VEC_OP_ASSIGN(<<=)
    _VEC_OP_ASSIGN(&=)
    _VEC_OP_ASSIGN(|=)
    _VEC_OP_ASSIGN(^=)

    // inline uint32_t operator[](size_t idx) const {
    //     return reinterpret_cast<const uint32_t *>(&v[0])[idx]; }
};


template<typename T>
struct GFi16v {
    T v;
    GFi16v& operator=(const T& val) { v = val; return *this; }
    inline GFi16v operator+(const GFi16v& other) const { return {add(v, other.v)}; }
    inline GFi16v operator+(uint32_t other) const { return {add(v, other)}; }
    inline GFi16v operator-(const GFi16v& other) const { return {sub(v, other.v)}; }
    inline GFi16v operator*(const GFi16v& other) const { return {mul(v, other.v)}; }
    inline GFi16v operator*(uint32_t other) const { return {mul(v, other)}; }
    inline GFi16v operator/(const GFi16v& other) const { return {div(v, other.v)}; }
    inline GFi16v operator-() const { return {neg(v)}; }
    inline uint32_t operator[](size_t idx) const { return v[idx]; }
    inline operator T const&() const { return v; }
    inline operator T&() { return v; }

    template<typename V>
    static inline T mul(T const& lhs, V const& rhs) {
        T res = lhs * rhs;
        T overflow = (res == 0 && lhs && rhs);
        res = (res & !overflow) + (1 & overflow);
        res = (res & 0xffff) - (res >> 16);
        res += (res >= 0x80000000) & 0x10001;
        return res;
    }

    template<typename V>
    static inline T add(T const& lhs, V const& rhs) {
        T res = lhs + rhs;
        res -= (res > 0x10001) & 0x10001;
        return res;
    }

    static inline T sub(T const& lhs, T const& rhs) {
        T res = lhs - rhs;
        res += (res >= 0x80000000) & 0x10001;
        return res;
    }

    template<typename V>
    static inline V neg(V const& a) {
        return 0x10001 - a;
    }
};


template<typename GFT>
class RSi16vImpl {
public:
    inline RSi16vImpl(GFi16 const& gf, size_t block_size, size_t ecc_len):
        gf(gf),
        block_size(block_size),
        ecc_len(ecc_len)
    {
        size_t nbits = __builtin_ctzl(block_size);

        uint32_t root = gf.exp(gf.div(gf.log(1), block_size));

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
            _rbo[i] = rbo16(i)  >> (16 - nbits);

        _ecc_mix_wv.resize(ecc_len);
        for (size_t j = 0; j < ecc_len; ++j) {
            auto w = _roots_iv[_rbo[block_size - ecc_len]][0];
            // w = - w ** (- rbo(i) * j) / block_size
            _ecc_mix_wv[j] = GFT{0} + gf.neg(gf.div(gf.pow(w, j), block_size));
        }
    }

    inline void encode(GFT block[]) {
        ct_butterfly(&_rootsv[0], block, block_size);

        for (size_t j = 0; j < ecc_len; ++j)
            block[j] = block[j] * _ecc_mix_wv[j];

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
                    block[start] = a + b;
                    block[start + stride] = a - b;
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // j = i - start
                    // GFT w = roots[exp_f * (i - start)];
                    GFT w = roots[exp_f * (i - start)];

                    // Cooley-Tukey butterfly
                    GFT a = block[i];
                    GFT b = block[i + stride];
                    GFT m = b * w;
                    block[i] = a + m;
                    block[i + stride] = a - m;
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
                    block[start] = a + b;
                    block[start + stride] = a - b;
                }
                for (size_t i = start + 1; i < start + stride; ++i) {
                    // Gentleman-Sande butterfly
                    // GFT w = roots[(i - start) << exp_f];
                    GFT w = roots[(i - start) << exp_f];

                    GFT a = block[i];
                    GFT b = block[i + stride];
                    block[i] = a + b;
                    block[i + stride] = (a - b) * w;
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
