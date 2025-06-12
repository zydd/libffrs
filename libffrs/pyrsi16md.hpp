/**************************************************************************
 * PyRSi16md.hpp
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

#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "util.hpp"

namespace py = pybind11;


typedef uint32_t u32x4 __attribute__((vector_size(4 * sizeof(uint32_t))));
typedef uint32_t u32x8 __attribute__((vector_size(8 * sizeof(uint32_t))));
typedef uint32_t u32x16 __attribute__((vector_size(16 * sizeof(uint32_t))));

template<typename T>
class GFi16v {
public:
    struct GFT {
        using GF = GFi16v;
        T v;
        GFT& operator=(const T& val) { v = val; return *this; }
        inline GFT operator+(const GFT& other) const { return {GFi16v::add(v, other.v)}; }
        inline GFT operator+(uint32_t other) const { return {GFi16v::add(v, other)}; }
        inline GFT operator-(const GFT& other) const { return {GFi16v::sub(v, other.v)}; }
        inline GFT operator*(const GFT& other) const { return {GFi16v::mul(v, other.v)}; }
        inline GFT operator*(uint32_t other) const { return {GFi16v::mul(v, other)}; }
        inline GFT operator/(const GFT& other) const { return {GFi16v::div(v, other.v)}; }
        inline GFT operator-() const { return {GFi16v::neg(v)}; }
        inline uint32_t operator[](size_t idx) const { return v[idx]; }
        inline operator T const&() const { return v; }
        inline operator T&() { return v; }
    };
    static_assert(sizeof(T) == sizeof(GFT));

    inline GFi16v(uint32_t primitive) {
        uint32_t x = 1;
        for (uint32_t i = 0; i < 0x10001; ++i) {
            _exp[i] = x;
            _log[x] = i;
            x = (x * primitive) % 0x10001;
        }
    }

    template<typename... Args>
    inline GFT operator()(Args... args) const {
        static_assert(sizeof...(Args) * sizeof(uint32_t) == sizeof(T));
        T v = {static_cast<uint32_t>(args)...};
        return {v};
    }

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

    inline uint32_t inv1(uint32_t a) const {
        return _exp[0x10000 - _log[a]];
    }

    inline uint32_t div1(uint32_t a, uint32_t b) const {
        if (a == 0)
            return 0;

        size_t r = _log[a] + 0x10000 - _log[b];
        if (r >= 0x10000)
            r -= 0x10000;

        return _exp[r];
    }

    inline uint32_t exp1(uint32_t a) const {
        return _exp[a];
    }

    inline uint32_t log1(uint32_t a) const {
        return _log[a];
    }

    inline uint32_t pow1(uint32_t a, uint32_t b) const {
        return _exp[(_log[a] * b) % 0x10000];
    }

private:
    std::array<uint32_t, 0x10001> _exp = {};
    std::array<uint32_t, 0x10001> _log = {};
};


template<typename GF>
class RSi16v {
public:
    using GFT = typename GF::GFT;

    inline RSi16v(size_t block_size, size_t ecc_len, uint32_t primitive):
        gf(primitive),
        block_size(block_size),
        ecc_len(ecc_len)
    {
        size_t nbits = __builtin_ctzl(block_size);

        uint32_t root = gf.exp1(gf.div1(gf.log1(1), block_size));

        if (root >= 0x8000)
            root = gf.neg(root);

        if (gf.pow1(root, block_size) != 1)
            throw std::runtime_error("Root of unity not found for block size");

        _rootsv.resize(block_size);
        _roots_iv.resize(block_size);
        uint32_t r = root;
        uint32_t r_i = gf.inv1(r);
        for (size_t i = 0; i < block_size; ++i) {
            // _roots[i] = r;
            // r = gf.pow(r, 2);
            _rootsv[i] = GFT{0} + gf.pow1(root, i);
            _roots_iv[i] = GFT{0} + gf.pow1(r_i, i);
        }

        _rbo.resize(block_size);
        for (size_t i = 0; i < block_size; ++i)
            _rbo[i] = rbo16(i)  >> (16 - nbits);

        _ecc_mix_wv.resize(ecc_len);
        for (size_t j = 0; j < ecc_len; ++j) {
            auto w = _roots_iv[_rbo[block_size - ecc_len]][0];
            // w = - w ** (- rbo(i) * j) / block_size
            _ecc_mix_wv[j] = GFT{0} + gf.neg(gf.div1(gf.pow1(w, j), block_size));
        }
    }

    inline void encode(GFT block[], size_t block_size) {
        ct_butterfly(&_rootsv[0], block, block_size);

        for (size_t j = 0; j < ecc_len; ++j)
            block[j] = block[j] * _ecc_mix_wv[j];

        for (size_t j = 1; j < block_size / ecc_len; ++j)
            std::copy_n(&block[0], ecc_len, &block[j * ecc_len]);
            // std::memcpy(&block[j * ecc_len], block, ecc_len * sizeof(GFT));

        gs_butterfly(ecc_len, &_roots_iv[0], block, block_size);
    }

protected:
    // GFT _nth_roots_of_unity[MaxFieldBits] = {0};
    GF gf;
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

using GFi16v4 = GFi16v<u32x4>;

class PyRSi16md : public RSi16v<GFi16v4> {
public:
    using GFT = GFi16v4::GFT;
    static constexpr size_t vec_size = sizeof(GFT) / sizeof(uint32_t);

    inline PyRSi16md(
            std::optional<size_t> block_size,
            std::optional<uint16_t> message_len,
            std::optional<uint16_t> ecc_len,
            uint32_t primitive):
        PyRSi16md(_get_ecc_len(block_size, message_len, ecc_len),
                // TODO: _get_block_size(block_size, message_len, ecc_len),
                block_size.value_or(256) / sizeof(uint16_t),
                primitive)
    {
        if (ecc_len && message_len) {
            if (block_size && *message_len + *ecc_len != *block_size) {
                throw py::value_error("block_size must be equal to message_len + ecc_len");
            } else if (!block_size) {
                set_block_size((size_t(*message_len) + size_t(*ecc_len)) / sizeof(uint16_t));
            }
        }
    }

    inline PyRSi16md(uint16_t ecc_len, size_t block_size, uint32_t primitive):
        RSi16v<GFi16v4>(block_size, ecc_len, primitive)
    {
        set_block_size(block_size);
    }

    inline void set_block_size(size_t block_size) {
        if (block_size <= ecc_len)
            throw py::value_error("block_size must be greater than ecc_len");

        if (block_size > 65536)
            throw py::value_error("block_size must be <= 65536");

        this->block_size = block_size;
    }

    template<typename T, typename U>
    inline void copy_blocks(const T *src, size_t src_size, U *dst) const {
        for (size_t j = 0; j < vec_size; ++j)
            for (size_t i = 0; i < src_size; ++i)
                dst[i + j * block_size] = src[i + j * src_size];
    }

    template<typename T, typename U>
    inline void copy_blocks_transposed(const T *src, size_t src_size, U *dst) const {
        for (size_t j = 0; j < vec_size; ++j)
            for (size_t i = 0; i < src_size; ++i)
                dst[i + j * block_size] = src[i * vec_size + j];
    }

    template<typename T, typename U>
    inline void copy_msg_transposed(const T *src, size_t src_size, U *dst) const {
        for (size_t j = 0; j < vec_size; ++j)
            for (size_t i = 0; i < src_size; ++i)
                dst[i * vec_size + j] = src[i + j * src_size];
    }

    inline void encode_block_v(const uint16_t src[], size_t msg_size, GFT temp[], uint16_t dst[]) {
        auto temp_u32 = reinterpret_cast<uint32_t *>(&temp[0]);
        copy_blocks(&src[0], msg_size, &dst[0]);

        copy_msg_transposed(&src[0], msg_size, &temp_u32[0]);
        std::fill_n(&temp[msg_size], ecc_len, GFT{0});
        // std::memset(&temp[msg_size], 0, ecc_len * sizeof(GFT));

        encode(&temp[0], block_size);
        copy_blocks_transposed(&temp_u32[0], ecc_len, &dst[msg_size]);
    }

    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        size_t msg_size = block_size - ecc_len;
        py_assert(buf.size == msg_size * vec_size);

        std::vector<GFT> temp(block_size);
        auto output = py::bytearray(nullptr, block_size * vec_size * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

        encode_block_v(&buf[0], msg_size, &temp[0], &output_data[0]);

        return output;
    }

    inline py::bytearray py_encode_blocks(buffer_ro<uint16_t> buf) {
        size_t msg_size = block_size - ecc_len;

        if (buf.size == 0 || block_size == 0 || block_size <= ecc_len)
            return {};

        size_t full_blocks = buf.size / msg_size;

        size_t output_size = full_blocks * (msg_size + ecc_len);

        // Last block will be smaller if input size is not divisible by msg_size
        size_t input_remainder = buf.size - full_blocks * msg_size;
        if (input_remainder > 0)
            output_size += input_remainder + ecc_len;

        auto output = py::bytearray(nullptr, output_size * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));
        std::vector<GFT> temp(block_size);
        
        size_t block = 0;
        if (full_blocks >= vec_size)
            for (; block < full_blocks; block += vec_size)
                encode_block_v(&buf[block * msg_size], msg_size, &temp[0], &output_data[block * block_size]);

        if (block != full_blocks)
            throw std::runtime_error("Number of blocks is not divisible by vector size: " + std::to_string(full_blocks));

        py_assert(input_remainder == 0, "Not implemented");
        if (input_remainder > 0) {
            auto output_block = &output_data[output_size - input_remainder - ecc_len];
            std::copy_n(&buf[buf.size - input_remainder], input_remainder, &output_block[0]);

            // std::fill_n(&temp[0], block_size, 0);
            // encode(&buf[buf.size - input_remainder], input_remainder, &temp[0], block_size);
            // std::copy_n(&temp[0], ecc_len, &output_block[input_remainder]);
        }
        return output;
    }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyRSi16md>(m, "RSi16md")
            .def_property_readonly("ecc_len", [](PyRSi16md& self) { return self.ecc_len * sizeof(uint16_t); })

            .def_property("block_size",
                [](PyRSi16md& self) { return self.block_size * sizeof(uint16_t); },
                &PyRSi16md::set_block_size)

            .def_property_readonly("message_len",
                [](PyRSi16md& self) { return (self.block_size - self.ecc_len) * sizeof(uint16_t); })

            // .def_property_readonly("gf", [](PyRSi16md& self) -> auto const& { return self.gf; })

            // .def_property_readonly("roots_of_unity", [](PyRSi16md& self) {
            //     return std::vector<GFT>(std::begin(self._nth_roots_of_unity), std::end(self._nth_roots_of_unity)); })

            .def(py::init<std::optional<uint16_t>, std::optional<uint16_t>, std::optional<uint16_t>, uint32_t>(), R"(
                Instantiate a Reed-Solomon encoder with the given configuration)",
                "block_size"_a = py::none(), "message_len"_a = py::none(), "ecc_len"_a = py::none(), "primitive"_a = 3)

            .def("__sizeof__", [](PyRSi16md& self) { return sizeof(self); })

            .def("encode", cast_args(&PyRSi16md::py_encode),
                R"(Systematic encode)",
                "buffer"_a)

            .def("encode_blocks", cast_args(&PyRSi16md::py_encode_blocks), R"(
                Encode the input buffer in blocks, storing the parity bytes right next to their corresponding block
                )",
                "buffer"_a)

            // .def("decode", cast_args(&PyRSi16md::py_decode),
            //     R"(Systematic decode)",
            //     "buffer"_a)

            .doc() = R"(Reed-Solomon coding over :math:`GF(2^8)`)";
    }

private:
    inline uint16_t _get_ecc_len(
            std::optional<size_t> block_size,
            std::optional<uint16_t> message_len,
            std::optional<uint16_t> ecc_len) {

        uint16_t res;
        if (ecc_len) {
            res = *ecc_len;
        } else if (message_len && block_size) {
            if (*message_len >= *block_size) {
                throw py::value_error("block_size must be greater than message_len");
            }
            res = *block_size - *message_len;
        } else {
            throw py::value_error("Must specify either (block_size, message_len) or ecc_len");
        }

        if (res % 2 != 0)
            throw py::value_error("ecc_len must be a multiple of 2: " + std::to_string(res));

        return res / sizeof(uint16_t);
    }
};
