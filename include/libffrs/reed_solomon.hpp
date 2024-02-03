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
#include <new>

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
class rs_encode_basic {
public:
    using GFT = typename GF::GFT;

    inline void encode(const GFT input[], size_t input_size, GFT output[]) {
        auto rs = RS::cast(this);

        // Multiply input by X^n, n = ecc_len+1  === append ecc_len bytes
        auto input_x_n = (GFT *) alloca(input_size + rs.ecc_len);
        std::copy_n(input, input_size, input_x_n);
        std::fill_n(&input_x_n[input_size], rs.ecc_len, 0x00);

        rs.gf.poly_mod(input_x_n, input_size + rs.ecc_len, rs.generator, rs.ecc_len + 1, output);

        if (rs.gf.prime != 2) {
            for (size_t i = 0; i < rs.ecc_len; ++i)
                output[i] = rs.gf.sub(0, output[i]);
        }
    }
};


template<typename GF, typename RS>
class rs_encode_basic_v2 {
public:
    using GFT = typename GF::GFT;

    inline void encode(const GFT input[], size_t input_size, GFT output[]) {
        auto rs = RS::cast(this);
        rs.gf.poly_mod_x_n(input, input_size, &rs.generator[1], rs.ecc_len, output);

        if (rs.gf.prime != 2) {
            for (size_t i = 0; i < rs.ecc_len; ++i)
                output[i] = rs.gf.sub(0, output[i]);
        }
    }
};


template<size_t MaxFieldElements, size_t MaxEccLen = MaxFieldElements>
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
                rs.gf.poly_mod(data, rs.ecc_len + 1, rs.generator, rs.ecc_len + 1, data);

                for (size_t j = 0; j < rs.ecc_len; ++j)
                    generator_lut[i][j] = data[j];
            }
        }

    private:
        alignas(size_t)
        uint8_t generator_lut[detail::align_size(MaxFieldElements, alignof(size_t))][detail::align_size(MaxEccLen, alignof(size_t))] = {};
    };
};


template<typename Word, size_t EccLen, size_t MaxFieldElements, size_t Stride>
struct rs_encode_slice_pw2 {
    static inline void encode(const void *lut, const uint8_t *input, size_t size, uint8_t *output) {
        size_t i = 0;
        WordB rem = {};
        auto ecc_len = EccLen;

        auto& generator_lut = *static_cast<const generator_lut_t *>(lut);

        if constexpr (Stride > 1) {
            for (; size - i >= Stride; i += Stride) {
                rem.word ^= *reinterpret_cast<const Word *>(&input[i]);

                Word t = 0;

                for (size_t j = 0; j < EccLen; ++j)
                    t ^= generator_lut[Stride-1 - j][rem.bytes[j]].word;

                for (size_t j = ecc_len; j < Stride; ++j)
                    t ^= generator_lut[Stride-1 - j][input[i + j]].word;

                rem.word = t;
            }
        }

        for (; i < size; ++i)
            rem.word = (rem.word >> 8) ^ generator_lut[0][rem.bytes[0] ^ input[i]].word;

        std::copy_n(rem.bytes, EccLen, output);
    }

    template<typename RS>
    static inline void init(RS& rs) {
        auto generator_lut = new generator_lut_t();

        for (size_t i = 0; i < 256; ++i) {
            generator_lut[0][i].bytes[0] = uint8_t(i);
            rs.gf.poly_mod_x_n(generator_lut[0][i].bytes, 1, &rs.generator[1], EccLen, generator_lut[0][i].bytes);
        }

        for (size_t i = 0; i < 256; ++i) {
            for (size_t j = 1; j < Stride; ++j) {
                Word prev = generator_lut[(j - 1)][i].word;
                generator_lut[j][i].word = (prev >> 8) ^ generator_lut[0][prev & 0xff].word;
            }
        }

        rs.generator_lut = generator_lut;
    }

private:
    union WordB { Word word; uint8_t bytes[sizeof(Word)]; };
    using generator_lut_t = WordB[Stride][MaxFieldElements];
};


template<size_t EccLen, size_t MaxFieldElements, size_t Stride, size_t Alignment = 32>
struct rs_encode_slice_generic_pw2 {
    static_assert(Stride >= EccLen);
    static_assert(Stride % EccLen == 0);

    static inline void encode(const void *lut, const uint8_t *input, size_t size, uint8_t *output) {
        size_t i = 0;
        Word rem = {};

        auto& generator_lut = *static_cast<const generator_lut_t *>(lut);

        if constexpr (Stride > 1) {
            for (; i + Stride <= size; i += Stride) {
                // rem.word ^= *reinterpret_cast<const Word *>(&input[i]);
                std::transform(&rem[0], &rem[EccLen], &input[i], &rem[0], std::bit_xor());

                Word t = {};

                for (size_t j = 0; j < EccLen; ++j)
                    // t ^= generator_lut[Stride-1 - j][rem.bytes[j]].word;
                    std::transform(&t[0], &t[EccLen], &generator_lut[Stride-1 - j][rem[j]][0], &t[0], std::bit_xor());

                for (size_t j = EccLen; j < Stride; ++j)
                    // t ^= generator_lut[Stride-1 - j][input[i + j]].word;
                    std::transform(&t[0], &t[EccLen], &generator_lut[Stride-1 - j][input[i + j]][0], &t[0], std::bit_xor());

                // rem.word = t;
                std::copy_n(&t[0], EccLen, &rem[0]);
            }
        }

        for (; i < size; ++i) {
            // rem.word = (rem.word >> 8) ^ generator_lut[0][rem.bytes[0] ^ input[i]].word;
            uint8_t rem0 = rem[0];
            std::rotate(&rem[0], &rem[1], &rem[EccLen]);
            rem[EccLen - 1] = 0;
            std::transform(&rem[0], &rem[EccLen], &generator_lut[0][rem0 ^ input[i]][0], &rem[0], std::bit_xor());
        }

        std::copy_n(&rem[0], EccLen, &output[0]);
    }

    template<typename RS>
    static inline void init(RS& rs) {
        auto generator_lut = new generator_lut_t();

        for (size_t i = 0; i < 256; ++i) {
            generator_lut[0][i][0] = uint8_t(i);
            rs.gf.poly_mod_x_n(&generator_lut[0][i][0], 1, &rs.generator[1], EccLen, &generator_lut[0][i][0]);
        }

        for (size_t i = 0; i < 256; ++i) {
            for (size_t j = 1; j < Stride; ++j) {
                // Word prev = generator_lut[j - 1][i].word;
                // generator_lut[j][i].word = (prev >> 8) ^ generator_lut[0][prev & 0xff].word;

                Word rem;
                std::copy_n(&generator_lut[j - 1][i][0], EccLen, &rem[0]);

                uint8_t rem0 = rem[0];
                std::rotate(&rem[0], &rem[1], &rem[EccLen]);
                rem[EccLen - 1] = 0;
                std::transform(&rem[0], &rem[EccLen], &generator_lut[0][rem0][0], &generator_lut[j][i][0], std::bit_xor());
            }
        }

        rs.generator_lut = generator_lut;
    }

private:
    struct alignas(Alignment) Word : std::array<uint8_t, detail::align_size(EccLen, Alignment)> { };
    using generator_lut_t = Word[Stride][MaxFieldElements];
};


template<size_t MaxEccLen>
struct rs_encode_slice_pw2_dispatch {
    template<typename GF, typename RS>
    class type {
    public:
        using GFT = typename GF::GFT;

        inline void encode(const uint8_t *input, size_t size, uint8_t *output) const {
            return _encode(generator_lut, input, size, output);
        }

    protected:
        void init() {
            auto& rs = RS::cast_mut(this);
            _dispatch.init[rs.ecc_len](rs);
            _encode = _dispatch.encode[rs.ecc_len];
        }

        template<typename, size_t, size_t, size_t>
        friend struct rs_encode_slice_pw2;
        template<size_t, size_t, size_t, size_t>
        friend struct rs_encode_slice_generic_pw2;

        void *generator_lut = {};

    private:
        using encode_fn_t = void(*)(const void *, const uint8_t *, size_t, uint8_t *);
        using init_fn_t = void(*)(RS&);
        encode_fn_t _encode = {};

        struct dispatch_t {
            constexpr dispatch_t() { fill(); }

            template<size_t EccLen = 1>
            constexpr void fill() {
                if constexpr(EccLen <= 2) {
                    encode[EccLen] = rs_encode_slice_pw2<uint32_t, EccLen, 256, 8>::encode;
                    init[EccLen] = rs_encode_slice_pw2<uint32_t, EccLen, 256, 8>::init;
                } else if constexpr(EccLen <= 8) {
                    encode[EccLen] = rs_encode_slice_pw2<uint64_t, EccLen, 256, 16>::encode;
                    init[EccLen] = rs_encode_slice_pw2<uint64_t, EccLen, 256, 16>::init;
#ifdef __GNUC__
                } else if constexpr(EccLen <= 16) {
                    encode[EccLen] = rs_encode_slice_pw2<__uint128_t, EccLen, 256, 16>::encode;
                    init[EccLen] = rs_encode_slice_pw2<__uint128_t, EccLen, 256, 16>::init;
#endif
                } else {
                    encode[EccLen] = rs_encode_slice_generic_pw2<EccLen, 256, EccLen>::encode;
                    init[EccLen] = rs_encode_slice_generic_pw2<EccLen, 256, EccLen>::init;
                }

                if constexpr(EccLen < MaxEccLen-1)
                    fill<EccLen + 1>();
            }

            encode_fn_t encode[MaxEccLen] = {};
            init_fn_t init[MaxEccLen] = {};
        };

        static constexpr dispatch_t _dispatch = {};
    };
};


template<size_t MaxFieldElements>
struct rs_synds_basic {
    template<typename GF, typename RS>
    class type {
    public:
        using synds_array_t = typename GF::GFT[MaxFieldElements];

        template<typename S, typename T>
        inline void synds(S const& data, size_t size, T const& rem, synds_array_t synds) const {
            auto rs = RS::cast(this);

            for (size_t i = 0; i < rs.ecc_len; ++i) {
                auto t = rs.gf.poly_eval(data, size, rs.generator_roots[i]);
                synds[i] = rs.gf.poly_eval(rem, rs.ecc_len, rs.generator_roots[i], t);
            }
        }
    };
};


template<typename Word, size_t MaxEccLen>
struct rs_synds_lut_pw2 {
    template<typename GF, typename RS>
    class type {
        static_assert(sizeof(typename GF::wide_mul_word_t) >= sizeof(Word));

    public:
        static constexpr size_t max_ecc_w = (MaxEccLen / sizeof(Word)) + !!(MaxEccLen % sizeof(Word));
        static constexpr size_t synds_size = max_ecc_w * sizeof(Word);
        using synds_array_t /* alignas(sizeof(Word)) */ = uint8_t[synds_size];
        using GFT = typename GF::GFT;

        inline void synds(const GFT *data, size_t size, const GFT *rem, synds_array_t synds) const {
            auto rs = RS::cast(this);
            auto synds_w = reinterpret_cast<Word *>(synds);
            for (size_t i = 0; i < ecc_w; ++i) {
                auto t = rs.gf.poly_eval_wide(data, size, generator_roots_wide[i]);
                synds_w[i] = Word(rs.gf.poly_eval_wide(rem, rs.ecc_len, generator_roots_wide[i], t));
            }
        }

    protected:
        inline void init() {
            auto rs = RS::cast(this);

            size_t i = 0;
            for (; i + sizeof(Word) <= rs.ecc_len; i += sizeof(Word)) {
                generator_roots_wide[i / sizeof(Word)] = reinterpret_cast<const Word *>(rs.generator_roots)[i / sizeof(Word)];
            }

            for (; i < rs.ecc_len; ++i) {
                reinterpret_cast<GFT *>(generator_roots_wide)[i] = rs.generator_roots[i];
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
        rs.synds(&data[0], size, rem, synds);

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
        rs.synds(&data[0], size, rem, synds);

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

    inline size_t berlekamp_massey(const GFT synds[/*ecc_len*/], GFT err_poly[/*ecc_len*/]) const {
        auto rs = RS::cast(this);
        auto prev = (GFT *) alloca(rs.ecc_len);
        std::fill_n(prev, rs.ecc_len, 0x00);
        auto temp = (GFT *) alloca(rs.ecc_len);
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
            const GFT synds[/*ecc_len*/], GFT err_poly[], const GFT err_pos[],
            const size_t err_count, GFT err_mag[]) const {
        auto rs = RS::cast(this);
        auto err_eval = (GFT *) alloca(rs.ecc_len * 2);

        auto synds_rev = (GFT *) alloca(rs.ecc_len);
        std::reverse_copy(synds, &synds[rs.ecc_len], synds_rev);

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
