/**************************************************************************
 * galois.hpp
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

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <functional>

#include "detail.hpp"

namespace ffrs {

using size_t = std::size_t;


template<typename T, template<class, class>typename...Fs>
class GF : public detail::CRTP<Fs<T, GF<T, Fs...>>...>  {
public:
    using detail::CRTP<Fs<T, GF<T, Fs...>>...>::CRTP;
};


template<typename T, typename GF>
class gf_data {
public:
    using GFT = T;
    const GFT prime;
    const GFT power;
    const GFT primitive;
    const GFT poly1;
    const size_t field_elements;

    inline gf_data(GFT prime, GFT power, GFT primitive, GFT poly1):
        prime(prime), power(power), primitive(primitive), poly1(poly1),
        field_elements(detail::ipow(prime, power))
    { }
};


/*************
 * Prime Field
 *************/


template<typename GFT, typename GF>
class gf_add_mod {
public:
    inline GFT add(GFT const& lhs, GFT const& rhs) const {
        auto& gf = GF::cast(this);
        // assert(gf.power == 1);
        return (lhs + rhs) % gf.prime;
    }

    inline GFT sub(GFT const& lhs, GFT const& rhs) const {
        auto& gf = GF::cast(this);
        // assert(gf.power == 1);
        if (lhs >= rhs)
            return (lhs - rhs) % gf.prime;
        else
            return ((gf.field_elements - rhs) + lhs) % gf.prime;
    }
};


template<typename GFT, typename GF>
class gf_mul_mod {
public:
    inline GFT mul(GFT const& lhs, GFT const& rhs) const {
        auto& gf = GF::cast(this);
        // assert(gf.power == 1);
        GFT res;
        if (__builtin_mul_overflow(lhs, rhs, &res))
            res = 1;
        return res % gf.prime;
    }
};


template<typename GFT, typename GF>
class gf_add_i16_shift {
public:
    inline GFT add(GFT const& lhs, GFT const& rhs) const {
        uint32_t res = lhs + rhs;
        if (res > 0x10001)
            res -= 0x10001;
        // res -= (res > 0x10001) * 0x10001;
        return res;
    }

    inline GFT sub(GFT const& lhs, GFT const& rhs) const {
        uint32_t res = lhs - rhs;
        if (int32_t(res) < 0)
            res += 0x10001;
        // res += (res >= 0x10000000) * 0x10001;
        return res;
    }

    inline GFT neg(GFT const& a) const {
        return 0x10001 - a;
    }
};


template<typename GFT, typename GF>
class gf_mul_i16_shift {
public:
    inline GFT mul(GFT const& lhs, GFT const& rhs) const {
        uint32_t res = lhs * rhs;
        if (res == 0 && lhs && rhs) [[unlikely]]
            return 1;

        // uint32_t overflow = (res == 0 && lhs && rhs) * 0xffffffff;
        // res = (res & ~overflow) + (1 & overflow);

        // uint32_t res;
        // if (__builtin_mul_overflow(lhs, rhs, &res))
        //     res = 1;

        res = (res & 0xffff) - (res >> 16);

        if (int32_t(res) < 0)
            res += 0x10001;

        return res;
    }
};


/*********
 * Generic
 *********/


template<template<class, class>typename MulOperation, size_t MaxFieldElements>
struct gf_mul_lut {
    template<typename GFT, typename GF>
    class type {
    public:
        static_assert(MaxFieldElements <= 4096); // limit 16MB

        inline GFT mul(GFT const& a, GFT const& b) const {
            auto& gf = GF::cast(this);
            return _mul[a * gf.field_elements + b];
        }

    protected:
        inline type() {
            auto& gf = GF::cast(this);
            using GF_mul = ffrs::GF<GFT, ffrs::gf_data, MulOperation>;
            auto gf_cpu = GF_mul(typename GF_mul::gf_data(gf.prime, gf.power, gf.primitive, gf.poly1));

            for (size_t i = 0; i < gf.field_elements; ++i) {
                for (size_t j = 0; j < gf.field_elements; ++j)
                    _mul[i * gf.field_elements + j] = gf_cpu.mul(i, j);
            }
        }
    private:
        GFT _mul[MaxFieldElements * MaxFieldElements] = {};
    };
};


template<template<class, class>typename MulOperation, size_t MaxFieldElements>
struct gf_exp_log_lut {
    template<typename GFT, typename GF>
    class type {
    public:
        inline GFT inv(GFT const& a) const {
            auto& gf = GF::cast(this);
            return _exp[gf.field_elements-1 - _log[a]];
        }

        inline GFT div(GFT const& a, GFT const& b) const {
            auto& gf = GF::cast(this);

            if (a == 0)
                return 0;

            size_t r = _log[a] + gf.field_elements-1 - _log[b];
            if (r >= gf.field_elements-1)
                r -= gf.field_elements-1;

            return _exp[r];
        }

        inline GFT exp(GFT const& a) const {
            return _exp[a];
        }

        inline GFT log(GFT const& a) const {
            return _log[a];
        }

        inline GFT pow(GFT const& a, GFT const& b) const {
            auto& gf = GF::cast(this);
            return _exp[(_log[a] * b) % (gf.field_elements - 1)];
        }

    protected:
        inline type() {
            auto& gf = GF::cast(this);
            using GF_mul = ffrs::GF<GFT, ffrs::gf_data, MulOperation>;
            auto gf_cpu = GF_mul(typename GF_mul::gf_data(gf.prime, gf.power, gf.primitive, gf.poly1));

            GFT x = 1;
            for (size_t i = 0; i < gf.field_elements; ++i) {
                _exp[i] = x;
                _log[x] = GFT(i);
                x = gf_cpu.mul(x, gf.primitive);
            }
        }

    private:
        std::array<GFT, MaxFieldElements> _exp = {};
        std::array<GFT, MaxFieldElements> _log = {};
    };
};


template<typename GFT, typename GF>
struct gf_poly {
    template<typename T, typename U>
    inline size_t ex_synth_div(T a[], size_t size_a, const U b[], size_t size_b) const {
        if (size_b > size_a)
            return 0;

        auto& gf = GF::cast(this);

        GFT norm_factor = gf.inv(b[0]);

        for (size_t i = 0; i < size_a - size_b + 1; ++i) {
            auto c = a[i] = gf.mul(a[i], norm_factor);

            if (c == 0)
                continue;

            for (size_t j = 1; j < size_b; ++j)
                a[i + j] = gf.sub(a[i + j], gf.mul(b[j], c));
        }

        return size_a - size_b + 1;
    }

    template<typename T>
    inline size_t poly_mod(
            const T a[], const size_t size_a,
            const T b[], const  size_t size_b,
            T rem[]) const {
        if (size_b < 2)
            return 0;

        const size_t size_r = size_b - 1;
        std::copy_n(a, std::min(size_a, size_r), rem);

        if (size_a < size_b)
            return size_a;

        auto& gf = GF::cast(this);

        GFT norm_factor = gf.inv(b[0]);
        b += 1;

        for (size_t i = 0; i < size_a - size_r; ++i) {
            T c = gf.mul(rem[0], norm_factor);
            std::rotate(rem, rem + 1, rem + size_r);
            rem[size_r - 1] = a[size_r + i];

            if (c == 0)
                continue;

            for (size_t j = 0; j < size_r; ++j)
                rem[j] = gf.sub(rem[j], gf.mul(b[j], c));
        }

        return size_r;
    }

    template<typename Ta, typename Tb, typename Tr>
    inline void poly_mod_x_n(
            const Ta a[], const size_t size_a,
            const Tb b[], const size_t size_b,
            Tr rem[]) const {
        auto& gf = GF::cast(this);

        if (size_a >= size_b) {
            std::copy_n(a, size_b, rem);
            for (size_t i = 0; i < size_a - size_b; ++i) {
                Tr c = rem[0];
                rem[0] = a[size_b + i];
                std::rotate(rem, rem + 1, rem + size_b);

                if (c == 0)
                    continue;

                for (size_t j = 0; j < size_b; ++j)
                    rem[j] = gf.sub(rem[j], gf.mul(b[j], c));
                // std::transform(&rem[0], &rem[size_b], &b[0], &rem[0],
                //         [c](T a, T b) { return gf.sub(a, gf.mul(b, c)); });
            }
        } else {
            std::copy_n(a, size_a, &rem[size_b - size_a]);
            std::fill_n(rem, size_b - size_a, 0x00);
        }

        for (size_t i = 0; i < size_b; ++i) {
            Tr c = rem[0];
            rem[0] = 0;
            std::rotate(rem, rem + 1, rem + size_b);

            if (c == 0)
                continue;

            for (size_t j = 0; j < size_b; ++j)
                rem[j] = gf.sub(rem[j], gf.mul(b[j], c));
        }
    }

    template<typename T>
    inline GFT poly_eval(const T poly[], const size_t size, GFT const& x, GFT r = 0) const {
        auto& gf = GF::cast(this);
        for (size_t i = 0; i < size; ++i)
            r = gf.add(gf.mul(r, x), poly[i]);
        return r;
    }

    template<typename T>
    inline void poly_shift(T poly[], const size_t size, const size_t n) const {
        size_t i = 0;
        for (; i < size - n; ++i)
            poly[i] = poly[i + n];
        for (; i < size; ++i)
            poly[i] = 0;
    }

    template<typename T>
    inline void poly_scale(T poly[], const size_t size, T const& a) const {
        auto& gf = GF::cast(this);
        for (size_t i = 0; i < size; ++i)
            poly[i] = gf.mul(a, poly[i]);
    }

    template<typename T>
    inline void poly_add(const T poly_a[], const T poly_b[], T output[], const size_t size) const {
        auto& gf = GF::cast(this);
        for (size_t i = 0; i < size; ++i)
            output[i] = gf.add(poly_a[i], poly_b[i]);
    }

    template<typename T>
    inline void poly_add(const T poly_a[], const size_t size_a, const T poly_b[], const size_t size_b, T output[]) const {
        auto& gf = GF::cast(this);

        if (size_a > size_b) {
            auto start = size_a - size_b;
            std::copy_n(poly_a, start, output);
            for (size_t i = start; i < size_a; ++i)
                output[i] = gf.add(poly_a[i], poly_b[i - start]);
        } else if (size_b > size_a) {
            auto start = size_b - size_a;
            std::copy_n(poly_b, start, output);
            for (size_t i = start; i < size_b; ++i)
                output[i] = gf.add(poly_a[i - start], poly_b[i]);
        } else {
            for (size_t i = 0; i < size_a; ++i)
                output[i] = gf.add(poly_a[i], poly_b[i]);
        }
    }

    template<typename T>
    inline void poly_sub(const T poly_a[], const T poly_b[], T output[], const size_t size) const {
        auto& gf = GF::cast(this);
        for (size_t i = 0; i < size; ++i)
            output[i] = gf.sub(poly_a[i], poly_b[i]);
    }

    template<typename T>
    inline void poly_sub(const T poly_a[], const size_t size_a, const T poly_b[], const size_t size_b, T output[]) const {
        auto& gf = GF::cast(this);

        if (size_a > size_b) {
            auto start = size_a - size_b;
            std::copy_n(poly_a, start, output);
            for (size_t i = start; i < size_a; ++i)
                output[i] = gf.sub(poly_a[i], poly_b[i - start]);
        } else if (size_b > size_a) {
            auto start = size_b - size_a;
            std::copy_n(poly_b, start, output);
            for (size_t i = start; i < size_b; ++i)
                output[i] = gf.sub(poly_a[i - start], poly_b[i]);
        } else {
            for (size_t i = 0; i < size_a; ++i)
                output[i] = gf.sub(poly_a[i], poly_b[i]);
        }
    }

    template<typename T, typename U>
    inline size_t poly_mul(
            const T poly_a[], const size_t size_a,
            const U poly_b[], const size_t size_b,
            GFT r[]) const {
        size_t res_len = size_a + size_b - 1;

        for (size_t i = 0; i < res_len; ++i)
            r[i] = 0;

        auto& gf = GF::cast(this);
        for (size_t i = 0; i < size_a; ++i)
            for (size_t j = 0; j < size_b; ++j)
                r[i + j] = gf.add(r[i + j], gf.mul(poly_a[i], poly_b[j]));

        return res_len;
    }
};


template<typename GFT, typename GF>
struct gf_poly_deriv_prime {
    template<typename T>
    inline size_t poly_deriv(T poly[], size_t size) const {
        auto& gf = GF::cast(this);
        for (size_t i = 1; i < size; ++i)
            poly[size - i] = gf.mul(poly[size - i - 1], i);
        poly[0] = 0;
        return size - 1;
    }
};


/********
 * Base 2
 ********/


template<typename GFT, typename GF>
struct gf_add_xor {
    inline GFT add(GFT const& lhs, GFT const& rhs) const { return lhs ^ rhs; }
    inline GFT sub(GFT const& lhs, GFT const& rhs) const { return lhs ^ rhs; }
};


template<typename GFT, typename GF>
class gf_mul_cpu_pw2 {
public:
    inline GFT mul(GFT const& a, GFT const& b) const {
        auto& gf = GF::cast(this);

        GFT r = 0;

        for (int i = iter; i >= 0; --i) {
            if (r & (gf.field_elements >> 1))
                r = (r << 1) ^ gf.poly1;
            else
                r = (r << 1);

            if (a & (1 << i))
                r ^= b;
        }

        return r;
    }

    inline gf_mul_cpu_pw2() {
        auto& gf = GF::cast(this);
        iter = int(detail::ilog2_floor(gf.field_elements >> 1));
    }
private:
    int iter;
};

template<typename GFT, typename GF>
class gf_mul_exp_log_lut {
public:
    inline GFT mul(GFT const& a, GFT const& b) const {
        if (a == 0 || b == 0)
            return 0;

        auto& gf = GF::cast(this);
        size_t r = gf.log(a) + gf.log(b);
        if (r >= gf.field_elements-1)
            r -= gf.field_elements-1;

        return gf.exp(GFT(r));
    }
};


template<typename Word>
struct gf_wide_mul {
    template<typename GFT, typename GF>
    class type {
        static_assert(std::is_same_v<GFT, uint8_t>);

    public:
        using wide_mul_word_t = Word;

        inline Word mul_wide(Word a, Word b) const {
            Word r = 0;
            for (int i = 7; i >= 0; --i) {
                Word m = r & rep_0x80;

                // static_assert((gf.poly1 & 0x80) == 0); // needed for the next statement
                m = m - (m >> 7);

                r = ((r & rep_0x7f) << 1) ^ (polyw & m);

                Word n = (a & (rep_0x01 << i)) >> i;
                n = (n << 8) - n;

                r ^= b & n;
            }
            return r;
        }

        inline Word poly_eval_wide(const uint8_t poly[], const size_t size, const Word x, Word r = 0) const {
            for (size_t i = 0; i < size; ++i)
                r = mul_wide(r, x) ^ (poly[i] * rep_0x01);
            return r;
        }

        static constexpr Word repeat_byte(uint8_t n) {
            Word p = 0;
            for (size_t i = 0; i < sizeof(Word); ++i)
                p |= Word(n) << (i * 8);
            return p;
        }

    protected:
        Word polyw;

        inline type() {
            auto& gf = GF::cast(this);
            assert((gf.poly1 & 0x80) == 0);
            polyw = repeat_byte(gf.poly1);
        }

        static constexpr Word rep_0x80 = repeat_byte(0x80);
        static constexpr Word rep_0x7f = repeat_byte(0x7f);
        static constexpr Word rep_0x01 = repeat_byte(0x01);
    };
};

template<typename GFT, typename GF>
struct gf_poly_deriv_pw2 {
    template<typename T>
    inline size_t poly_deriv(T poly[], size_t size) const {
        for (size_t i = 1; i < size; ++i)
            poly[size - i] = (i & 1) ? poly[size - i - 1] : 0;

        return size - 1;
    }
};


}
