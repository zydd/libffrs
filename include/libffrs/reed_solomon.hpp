/**************************************************************************
 * reed_solomon.hpp
 *
 * Copyright 2024 Gabriel Machado
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
#include <functional>

#include "detail.hpp"
#include "galois.hpp"

namespace ffrs {


template<typename GF, template<class, class>typename...Fs>
class RS : public detail::CRTP<RS<GF, Fs...>, Fs<GF, RS<GF, Fs...>>...>  {
public:
    using GFT = typename GF::GFT;
    const GF gf;
    const GFT ecc_len;

    inline RS(GF const& gf, GFT ecc_len):
        gf(gf), ecc_len(ecc_len)
    {
        detail::CRTP<RS<GF, Fs...>, Fs<GF, RS<GF, Fs...>>...>::init();
    }
};


template<size_t MaxEccLen>
struct rs_generator {
    template<typename GF, typename RS>
    class type {
    public:
        using GFT = typename GF::GFT;
        GFT generator[MaxEccLen + 1] = {};
        GFT generator_roots[MaxEccLen] = {};

    protected:
        inline void init() {
            auto rs = RS::cast(this);

            auto temp = (GFT *) alloca(rs.ecc_len + 1);
            std::fill_n(temp, rs.ecc_len + 1, 0x00);

            auto p1 = (rs.ecc_len & 1) ? &generator[0] : &temp[0];
            auto p2 = (rs.ecc_len & 1) ? &temp[0] : &generator[0];

            size_t len = 1;
            p2[0] = 1;

            for (GFT i = 0; i < rs.ecc_len; ++i) {
                generator_roots[i] = rs.gf.exp(i);
                GFT factor[] = {1, rs.gf.sub(0, generator_roots[i])};
                len = rs.gf.poly_mul(p2, len, factor, 2, p1);

                auto t = p1;
                p1 = p2;
                p2 = t;
            }
        }
    };
};

template<typename GF, typename RS>
struct rs_encode_basic {
    template<typename In, typename Out>
    inline void encode(In const& input, size_t input_size, Out output) {
        auto rs = RS::cast(this);

        // Multiply input by X^n, n = ecc_len+1  === append ecc_len bytes
        std::remove_reference_t<decltype(input[0])> input_x_n[input_size + rs.ecc_len] = {};
        std::copy_n(input, input_size, input_x_n);

        rs.gf.poly_mod(input_x_n, input_size + rs.ecc_len, rs.generator, rs.ecc_len + 1, output);

        if (rs.gf.prime != 2) {
            for (size_t i = 0; i < rs.ecc_len; ++i)
                output[i] = rs.gf.sub(0, output[i]);
        }
    }
};


template<typename GF, typename RS>
struct rs_encode_basic_v2 {
    template<typename In, typename Out>
    inline void encode(In const& input, size_t input_size, Out output) {
        auto rs = RS::cast(this);
        rs.gf.poly_mod_x_n(input, input_size, &rs.generator[1], rs.ecc_len, output);

        if (rs.gf.prime != 2) {
            for (size_t i = 0; i < rs.ecc_len; ++i)
                output[i] = rs.gf.sub(0, output[i]);
        }
    }
};


template<size_t MaxFieldElements>
struct rs_encode_lut_pw2 {
    template<typename GF, typename RS>
    class type {
    public:
        inline void encode(const uint8_t *input, size_t size, uint8_t *output) const {
            auto rs = RS::cast(this);
            std::fill_n(output, rs.ecc_len, 0x00);
            for (size_t i = 0; i < size; ++i) {
                uint8_t pos = output[0] ^ input[i];
                output[0] = 0;
                std::rotate(output, output + 1, output + rs.ecc_len);
                std::transform(output, output + rs.ecc_len,
                        &generator_lut[pos][0],
                        output, std::bit_xor());
            }
        }

    protected:
        inline void init() {
            auto rs = RS::cast(this);
            for (size_t i = 0; i < rs.gf.field_elements; ++i) {
                auto data = (uint8_t *) alloca(rs.ecc_len + 1);
                std::fill_n(data, rs.ecc_len + 1, 0x00);
                data[0] = uint8_t(i);
                rs.gf.ex_synth_div(&data[0], rs.ecc_len + 1, &rs.generator[0], rs.ecc_len + 1);

                for (size_t j = 0; j < rs.ecc_len; ++j)
                    generator_lut[i][j] = data[j + 1];
            }
        }

    private:
        uint8_t generator_lut[MaxFieldElements][MaxFieldElements] = {};
    };
};


template<typename Word, size_t MaxFieldElements, size_t Stride = sizeof(Word)>
struct rs_encode_slice_pw2 {
    template<typename GF, typename RS>
    class type {
        public:
        using GFT = typename GF::GFT;

        inline void encode(const uint8_t *input, size_t size, uint8_t *output) {
            auto rs = RS::cast(this);
            size_t i = 0;
            Word rem = 0;

            if constexpr (Stride > 1) {
                for (; size - i >= Stride; i += Stride) {
                    auto in = reinterpret_cast<const Word *>(&input[i]);
                    rem ^= *in++;
                    Word t = 0;
                    for (size_t j = 0; j < rs.ecc_len; ++j)
                        t ^= generator_lut[Stride - j - 1][(rem >> (8 * j)) & 0xff];

                    for (size_t k = 1; k < Stride / rs.ecc_len; ++k) {
                        auto two = *in++;
                        for (size_t j = 0; j < rs.ecc_len; ++j)
                            t ^= generator_lut[Stride - k * rs.ecc_len - j - 1][(two >> (8 * j)) & 0xff];
                    }
                    rem = t;
                }
            }

            for (; i < size; ++i)
                rem = (rem >> 8) ^ generator_lut[0][(rem & 0xff) ^ input[i]];

            // endianess dependent
            for (size_t i = 0; i < rs.ecc_len; ++i)
                output[i] = uint8_t(rem >> (8 * i));
        }

    protected:
        inline void init() {
            auto rs = RS::cast(this);

            assert(rs.ecc_len == sizeof(Word));
            assert(Stride == 1 || Stride % rs.ecc_len == 0);

            for (size_t i = 0; i < rs.gf.field_elements; ++i) {
                uint8_t data[sizeof(Word) + 1] = {0};
                data[0] = uint8_t(i);
                rs.gf.ex_synth_div(&data[0], sizeof(Word) + 1, &rs.generator[0], sizeof(Word) + 1);

                // endianess dependent
                for (size_t j = 0; j < sizeof(Word); ++j)
                    generator_lut[0][i] |= Word(data[j + 1]) << j * 8;
            }

            for (size_t i = 0; i < rs.gf.field_elements; ++i) {
                for (size_t j = 1; j < Stride; ++j)
                    generator_lut[j][i] =
                            (generator_lut[j - 1][i] >> 8) ^
                            generator_lut[0][generator_lut[j - 1][i] % rs.gf.field_elements];
            }
        }

    private:
        Word generator_lut[Stride][MaxFieldElements] = {};
    };
};


template<size_t MaxFieldElements>
struct rs_synds_basic {
    template<typename GF, typename RS>
    class type {
    public:
        using synds_array_t = typename GF::GFT[MaxFieldElements];

        template<typename S, typename T>
        inline void synds(synds_array_t synds, S const& data, size_t size, T const& rem) const {
            auto rs = RS::cast(this);
            for (size_t i = 0; i < rs.ecc_len; ++i) {
                auto t = rs.gf.poly_eval(data, size, rs.generator_roots[i]);
                synds[rs.ecc_len - i - 1] = rs.gf.poly_eval(rem, rs.ecc_len, rs.generator_roots[i], t);
            }
        }
    };
};


template<typename Word, size_t MaxEccLen>
struct rs_synds_lut_pw2 {
    template<typename GF, typename RS>
    class type {
    public:
        static constexpr size_t max_ecc_w = (MaxEccLen / sizeof(Word)) + !!(MaxEccLen % sizeof(Word));
        static constexpr size_t synds_size = max_ecc_w * sizeof(Word);
        using synds_array_t /* alignas(sizeof(Word)) */ = uint8_t[synds_size];

        inline void synds(uint8_t synds[], const uint8_t *data, size_t size, const uint8_t *rem) {
            auto rs = RS::cast(this);
            auto synds_w = reinterpret_cast<Word *>(synds);
            for (size_t i = 0; i < ecc_w; ++i) {
                auto t = rs.gf.poly_eval_wide(data, size, generator_roots_wide[i]);

                static_assert(sizeof(typename GF::wide_mul_word_t) >= sizeof(Word));

                synds_w[i] = Word(rs.gf.poly_eval_wide(rem, rs.ecc_len, generator_roots_wide[i], t));
            }
        }

    protected:
        inline void init() {
            auto rs = RS::cast(this);
            for (auto i = 0; i < rs.ecc_len; ++i) {
                // reverse-order syndromes, endianess dependent
                generator_roots_wide[i / sizeof(Word)] |= Word(rs.gf.exp(rs.ecc_len - i - 1)) << (i % sizeof(Word)) * 8;
            }

            ecc_w = (rs.ecc_len / sizeof(Word)) + !!(rs.ecc_len % sizeof(Word));
        }

    private:
        Word generator_roots_wide[max_ecc_w] = {};
        size_t ecc_w = {};
    };
};


template<typename GF, typename RS>
class rs_roots_eval_basic {
public:
    using GFT = typename GF::GFT;

    inline size_t roots(
            const GFT poly[], const size_t poly_size,
            const size_t max_search_pos,
            GFT roots[]) const {
        auto rs = RS::cast(this);
        size_t count = 0;

        for (size_t i = 0; i < max_search_pos; ++i) {
            if (rs.gf.poly_eval(poly, poly_size, rs.gf.inv(rs.gf.exp(i))) == 0)
                roots[count++] = i;
        }

        return count;
    }
};


template<typename GF, typename RS>
class rs_roots_eval_uint8_chien {
public:
    inline size_t roots(
            const uint8_t poly[], const size_t poly_size,
            const size_t /*max_search_pos*/,
            uint8_t roots[]) {
        auto rs = RS::cast(this);
        size_t count = 0;

        auto coefs = (uint8_t *) alloca(rs.ecc_len);
        std::reverse_copy(&poly[0], &poly[poly_size], &coefs[0]);

        for (int i = 254; i >= 0; --i) {
            uint8_t sum = 1;

            for (size_t j = 1; j < poly_size; ++j) {
                coefs[j] = rs.gf.mul(coefs[j], rs.gf.exp(j));
                sum = rs.gf.add(sum, coefs[j]);
            }

            if (sum == 0) {
                roots[count] = i;
                if (++count >= poly_size-1)
                    break;
            }
        }

        return count;
    }
};


template<typename Word>
struct rs_roots_eval_lut_pw2 {
    template<typename GF, typename RS>
    class type {
    public:
        inline size_t roots(
                const uint8_t poly[], const size_t poly_size,
                const size_t max_search_pos,
                uint8_t roots[]) {
            auto rs = RS::cast(this);
            size_t count = 0;

            for (size_t i = 0; i <= max_search_pos/sizeof(Word); ++i) {
                Word eval = rs.gf.poly_eval_wide(poly, poly_size, err_poly_roots.word[i]);
                for (size_t j = 0; j < sizeof(Word); ++j) {
                    if (reinterpret_cast<uint8_t *>(&eval)[j] == 0) {
                        size_t pos = i * sizeof(Word) + j;

                        if (pos < max_search_pos)
                            roots[count++] = uint8_t(pos);
                        else
                            break;
                    }
                }
            }

            return count;
        }
    protected:
        inline void init() {
            auto rs = RS::cast(this);
            for (unsigned i = 0; i < 256; ++i)
                err_poly_roots.u8[i] = rs.gf.inv(rs.gf.exp(uint8_t(i)));
        }

    private:
        union {
            std::array<uint8_t, 256> u8;
            std::array<Word, 256 / sizeof(Word)> word;
        } err_poly_roots = {};
    };
};


template<typename GF, typename RS>
struct rs_decode {
    using GFT = typename GF::GFT;

    template<typename T, typename U>
    inline bool decode(T data, const size_t size, U rem) const {
        auto rs = RS::cast(this);
        typename RS::synds_array_t synds;
        rs.synds(synds, &data[0], size, rem);

        if (std::all_of(&synds[0], &synds[rs.ecc_len], std::logical_not()))
            return true;

        auto err_poly = (GFT *) alloca(rs.ecc_len);
        auto errors = berlekamp_massey(synds, err_poly);

        if (errors > rs.ecc_len / 2)
            return false;

        auto err_pos = (GFT *) alloca(rs.ecc_len / 2);
        auto roots = rs.roots(&err_poly[rs.ecc_len-errors-1], errors+1, size + rs.ecc_len, err_pos);

        if (errors != roots)
            return false;

        auto err_mag = (GFT *) alloca(rs.ecc_len / 2);
        forney(synds, &err_poly[rs.ecc_len-errors-1], err_pos, errors, err_mag);

        for (size_t i = 0; i < errors; ++i) {
            size_t pos = size + rs.ecc_len - 1 - err_pos[i];
            if (pos >= size + rs.ecc_len)
                return false;

            if (pos < size)
                data[pos] = rs.gf.add(data[pos], err_mag[i]);
            else
                rem[pos - size] = rs.gf.add(rem[pos - size], err_mag[i]);
        }

        return true;
    }

    template<typename T, typename U, typename V>
    inline bool decode(T data, size_t size, U rem, const V err_idx, size_t errors) const {
        auto rs = RS::cast(this);
        if (errors > rs.ecc_len)
            return false;

        typename RS::synds_array_t synds;
        rs.synds(synds, &data[0], size, rem);

        if (std::all_of(&synds[0], &synds[rs.ecc_len], std::logical_not()))
            return true;

        auto err_pos = (GFT *) alloca(rs.ecc_len);
        for (size_t i = 0; i < errors; ++i) {
            if (err_idx[i] > size + rs.ecc_len - 1)
                return false;

            err_pos[i] = size + rs.ecc_len - 1 - err_idx[i];
        }

        auto err_poly = (GFT *) alloca(rs.ecc_len + 1);
        err_poly[0] = 1;
        size_t err_poly_len = 1;

        {
            auto temp = (GFT *) alloca(rs.ecc_len + 1);
            std::fill_n(temp, rs.ecc_len + 1, 0x00);
            temp[0] = 1;

            auto p1 = (errors & 1) ? &err_poly[0] : &temp[0];
            auto p2 = (errors & 1) ? &temp[0] : &err_poly[0];

            for (size_t i = 0; i < errors; ++i) {
                GFT factor[2] = {rs.gf.sub(0, rs.gf.exp(err_pos[i])), 1};
                err_poly_len = rs.gf.poly_mul(p2, err_poly_len, factor, 2, p1);
                std::swap(p1, p2);
            }
        }

        assert(err_poly_len == errors + 1);

        auto err_mag = (GFT *) alloca(rs.ecc_len);
        forney(synds, err_poly, err_pos, errors, err_mag);

        for (size_t i = 0; i < errors; ++i) {
            size_t pos = err_idx[i];
            if (pos >= size + rs.ecc_len)
                return false;

            if (pos < size)
                data[pos] = rs.gf.add(data[pos], err_mag[i]);
            else
                rem[pos - size] = rs.gf.add(rem[pos - size], err_mag[i]);
        }

        return true;
    }

    inline size_t berlekamp_massey(const GFT synds_rev[/*ecc_len*/], GFT err_poly[/*ecc_len*/]) const {
        auto rs = RS::cast(this);
        auto prev = (GFT *) alloca(rs.ecc_len);
        std::fill_n(prev, rs.ecc_len, 0x00);
        auto temp = (GFT *) alloca(rs.ecc_len);
        auto synds = (GFT *) alloca(rs.ecc_len);
        std::reverse_copy(synds_rev, &synds_rev[rs.ecc_len], synds);
        std::fill_n(err_poly, rs.ecc_len, 0);

        prev[rs.ecc_len-1] = 1;
        err_poly[rs.ecc_len-1] = 1;

        size_t errors = 0;
        size_t m = 1;
        GFT b = 1;

        for (size_t n = 0; n < rs.ecc_len; ++n) {
            GFT d = synds[n]; // discrepancy
            for (size_t i = 1; i < errors + 1; ++i)
                d = rs.gf.add(d, rs.gf.mul(err_poly[rs.ecc_len - 1 - i], synds[n-i]));

            if (d == 0) {  // discrepancy is already zero
                m = m + 1;
            } else if (2 * errors <= n) {
                std::copy_n(err_poly, rs.ecc_len, temp);

                rs.gf.poly_shift(prev, rs.ecc_len, m);
                rs.gf.poly_scale(prev, rs.ecc_len, rs.gf.div(d, b));

                rs.gf.poly_sub(err_poly, prev, err_poly, rs.ecc_len);

                errors = n + 1 - errors;
                std::copy_n(temp, rs.ecc_len, prev);

                b = d;
                m = 1;
            } else {
                std::copy_n(prev, rs.ecc_len, temp);

                rs.gf.poly_shift(temp, rs.ecc_len, m);

                rs.gf.poly_scale(temp, rs.ecc_len, rs.gf.div(d, b));
                rs.gf.poly_sub(err_poly, temp, err_poly, rs.ecc_len);

                m = m + 1;
            }
        }

        return errors;
    }

    inline void forney(
            const GFT synds_rev[/*ecc_len*/], GFT err_poly[], const GFT err_pos[],
            const size_t err_count, GFT err_mag[]) const {
        auto rs = RS::cast(this);
        auto err_eval = (GFT *) alloca(rs.ecc_len * 2);
        auto err_eval_size = rs.gf.poly_mul(
                synds_rev, rs.ecc_len,
                err_poly, err_count + 1,
                err_eval);

        auto x_poly = (GFT *) alloca(rs.ecc_len + 1);
        std::fill_n(x_poly, rs.ecc_len + 1, 0x00);
        x_poly[0] = 1;
        auto err_eval_begin = rs.gf.ex_synth_div(
                err_eval, err_eval_size,
                x_poly, rs.ecc_len + 1);
        while (err_eval[err_eval_begin] == 0) {
            err_eval_begin++;
            assert(err_eval_begin < err_eval_size);
        }
        err_eval_size = err_eval_size - err_eval_begin;

        rs.gf.poly_deriv(err_poly, err_count + 1);
        auto err_poly_deriv = err_poly + 1;
        auto err_poly_deriv_size = err_count;

        for (size_t i = 0; i < err_count; ++i) {
            auto xi = rs.gf.exp(err_pos[i]);
            auto xi_inv = rs.gf.inv(xi);

            auto n = rs.gf.poly_eval(&err_eval[err_eval_begin], err_eval_size, xi_inv);
            auto d = rs.gf.poly_eval(err_poly_deriv, err_poly_deriv_size, xi_inv);

            err_mag[i] = rs.gf.mul(xi, rs.gf.div(n, d));
        }
    }
};

}
