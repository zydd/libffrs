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

#include "galois.hpp"

namespace ffrs {


template<typename GF, template<class, class>typename...Fs>
class NTT : public detail::CRTP<Fs<GF, NTT<GF, Fs...>>...>  {
public:
    using detail::CRTP<Fs<GF, NTT>...>::CRTP;
};

template<typename GF, typename RS>
class ntt_data {
public:
    using GFT = typename GF::GFT;
    const GF gf;
    const GFT root_of_unity;

    inline ntt_data(GF&& gf, GFT root_of_unity):
        gf(std::move(gf)), root_of_unity(root_of_unity)
    { }
};


template<typename Word>
struct ntt_eval {
    template<typename GF, typename NTT>
    class type {
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

        inline Word ntt(const uint8_t data[], const size_t size, Word r = 0) const {
            for (size_t i = 0; i < size; ++i)
                r = mul_wide(r, roots) ^ (data[i] * rep_0x01);
            return r;
        }

        inline type() {
            auto& ntt = NTT::cast(this);
            assert((ntt.gf.poly1 & 0x80) == 0);
            polyw = _repeat_byte(ntt.gf.poly1);

            roots = 0;
            for (size_t i = 0; i < sizeof(Word); ++i)
                roots |= Word(ntt.gf.pow(ntt.root_of_unity, i)) << (i * 8);
        }

    private:
        Word polyw;
        Word roots;

        static constexpr Word _repeat_byte(uint8_t n) {
            Word p = 0;
            for (size_t i = 0; i < sizeof(Word); ++i)
                p |= Word(n) << (i * 8);
            return p;
        }
        static constexpr Word rep_0x80 = _repeat_byte(0x80);
        static constexpr Word rep_0x7f = _repeat_byte(0x7f);
        static constexpr Word rep_0x01 = _repeat_byte(0x01);
    };
};


}
