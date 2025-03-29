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

template<typename GF, typename NTT>
class ntt_data {
public:
    using GFT = typename GF::GFT;
    const GF gf;
    const GFT root_of_unity;

    inline ntt_data(GF&& gf, GFT root_of_unity):
        gf(std::move(gf)), root_of_unity(root_of_unity)
    { }
};


template<typename GF, typename NTT>
class ntt_eval {
public:
    using Word = typename GF::wide_mul_word_t;

    inline Word ntt(const uint8_t data[], const size_t size, Word r = 0) const {
        auto& ntt = NTT::cast(this);
        for (size_t i = 0; i < size; ++i)
            r = ntt.gf.mul_wide(r, _ntt_eval_roots) ^ (data[i] * rep_0x01);
        return r;
    }

    inline ntt_eval() {
        auto& ntt = NTT::cast(this);
        assert((ntt.gf.poly1 & 0x80) == 0);
        _ntt_eval_polyw = GF::repeat_byte(ntt.gf.poly1);

        _ntt_eval_roots = 0;
        for (size_t i = 0; i < sizeof(Word); ++i)
            _ntt_eval_roots |= Word(ntt.gf.pow(ntt.root_of_unity, i)) << (i * 8);
    }

private:
    Word _ntt_eval_polyw;
    Word _ntt_eval_roots;
    Word rep_0x01 = GF::repeat_byte(0x01);
};


template<size_t MaxFieldElements>
struct ntt_eval_lut {
    template<typename GF, typename NTT>
    class type {
    public:
        using GFT = typename GF::GFT;
        using Word = typename GF::wide_mul_word_t;

        inline Word ntt(const uint8_t data[], const size_t size, Word r = 0) const {
            auto& ntt = NTT::cast(this);
            for (size_t i = 0; i < size; ++i)
                r ^= ntt._ntt_eval_term[i][data[i]];
            return r;
        }

        inline type() {
            auto& ntt = NTT::cast(this);
            for (size_t x = 0; x < MaxFieldElements; ++x) {
                for (size_t i = 0; i < MaxFieldElements; ++i) {
                    _ntt_eval_term[i][x] = 0;
                    for (size_t j = 0; j < sizeof(Word) / sizeof(GFT); ++j) {
                        GFT a = ntt.gf.pow(ntt.root_of_unity, i * j);
                        // reinterpret_cast<GFT *>(_ntt_eval_term[i])[j] = ntt.gf.mul(a, i);
                        _ntt_eval_term[i][x] |= Word(ntt.gf.mul(a, x)) << (j * 8);
                    }
                }
            }
        };

    private:
        Word _ntt_eval_term[MaxFieldElements][MaxFieldElements];
    };
};


}
