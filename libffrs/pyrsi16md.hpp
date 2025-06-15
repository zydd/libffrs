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
#include "rsi16md_impl.hpp"
#include "rsi16md.h"

namespace py = pybind11;


class PyRSi16md {
private:
    struct rs_data {
        size_t block_size;
        size_t ecc_len;
    };

    GFi16 gf;
    RSi16vImpl<uint32_t> rs16;
    RSi16v<4> rs16v4;
    RSi16v<8> rs16v8;
    RSi16v<16> rs16v16;
    const size_t block_size;
    const size_t msg_size;
    const size_t ecc_len;
    bool simd_x16;
    bool simd_x8;
    bool simd_x4;
    bool inline_ecc;

    inline PyRSi16md(
        rs_data&& args,
        uint32_t primitive,
        bool inline_ecc,
        std::optional<bool> simd_x16,
        std::optional<bool> simd_x8,
        std::optional<bool> simd_x4
    ):
        gf(GFi16::gf_data(0x10001, 1, primitive, 0)),
        rs16(gf, args.block_size, args.ecc_len),
        rs16v4(gf, args.block_size, args.ecc_len),
        rs16v8(gf, args.block_size, args.ecc_len),
        rs16v16(gf, args.block_size, args.ecc_len),
        block_size(args.block_size),
        msg_size(args.block_size - args.ecc_len),
        ecc_len(args.ecc_len),
        simd_x16(simd_x16.value_or(__builtin_cpu_supports("avx512f"))),
        simd_x8(simd_x8.value_or(__builtin_cpu_supports("avx2"))),
        simd_x4(simd_x4.value_or(__builtin_cpu_supports("sse2"))),
        inline_ecc(inline_ecc)
    { }

public:
    static constexpr size_t max_vec_size = 16;
    static constexpr size_t vec_align = max_vec_size * sizeof(uint32_t);

    inline PyRSi16md(
            std::optional<size_t> block_size,
            std::optional<uint16_t> message_len,
            std::optional<uint16_t> ecc_len,
            uint32_t primitive,
            bool inline_ecc,
            std::optional<bool> simd_x16,
            std::optional<bool> simd_x8,
            std::optional<bool> simd_x4
    ):
        PyRSi16md(
            _get_args(block_size, message_len, ecc_len),
            primitive,
            inline_ecc,
            simd_x16,
            simd_x8,
            simd_x4
        )
    { }

    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        constexpr size_t vec_size = 1;
        py_assert(buf.size == msg_size * vec_size);

        size_t output_size = (inline_ecc ? msg_size + ecc_len : ecc_len);

        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_size * vec_size]);
        auto output = py::bytearray(nullptr, output_size * vec_size * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

        if (inline_ecc)
            encode_block_inline<vec_size>(rs16, &buf[0], &temp[0], &output_data[0]);
        else
            encode_block_external<vec_size>(rs16, &buf[0], &temp[0], &output_data[0]);

        return output;
    }

    inline py::bytearray py_encode_blocks(buffer_ro<uint16_t> buf) {
        if (buf.size == 0 || block_size == 0 || block_size <= ecc_len)
            return {};

        size_t full_blocks = buf.size / msg_size;
        size_t output_size = full_blocks * (inline_ecc ? msg_size + ecc_len : ecc_len);

        // Last block will be smaller if input size is not divisible by msg_size
        size_t input_remainder = buf.size - full_blocks * msg_size;
        if (input_remainder > 0)
            output_size += inline_ecc ? input_remainder + ecc_len : ecc_len;

        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_size * max_vec_size]);
        auto output = py::bytearray(nullptr, output_size * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

        size_t block;
        if (inline_ecc)
            block = encode_blocks_inline(&buf[0], full_blocks, &temp[0], &output_data[0]);
        else
            block = encode_blocks_external(&buf[0], full_blocks, &temp[0], &output_data[0]);

        if (block != full_blocks)
            throw std::runtime_error("Number of blocks is not divisible by vector size: " +  std::to_string(block) + "/" + std::to_string(full_blocks));

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

            .def_property_readonly("block_size",
                [](PyRSi16md& self) { return self.block_size * sizeof(uint16_t); })

            .def_property_readonly("message_len",
                [](PyRSi16md& self) { return (self.block_size - self.ecc_len) * sizeof(uint16_t); })

            .def(py::init<
                    std::optional<uint16_t>,
                    std::optional<uint16_t>,
                    std::optional<uint16_t>,
                    uint32_t,
                    bool,
                    std::optional<bool>,
                    std::optional<bool>,
                    std::optional<bool>>(),
                R"(Instantiate a Reed-Solomon encoder with the given configuration)",
                "block_size"_a = py::none(),
                "message_len"_a = py::none(),
                "ecc_len"_a = py::none(),
                "primitive"_a = 3,
                "inline"_a = false,
                "simd_x16"_a = py::none(),
                "simd_x8"_a = py::none(),
                "simd_x4"_a = py::none()
            )

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
    inline rs_data _get_args(
            std::optional<size_t> block_size,
            std::optional<uint16_t> message_len,
            std::optional<uint16_t> ecc_len) {

        size_t block, message, ecc;
        if (ecc_len && block_size) {
            ecc = *ecc_len;
            block = *block_size;
            message = block - ecc;

            if (message_len and *message_len != block - ecc)
                throw py::value_error("block_size = message_len + ecc_len");
        } else if (message_len && block_size) {
            if (*message_len >= *block_size)
                throw py::value_error("block_size must be greater than message_len");

            block = *block_size;
            message = *message_len;
            ecc = block - message;
        } else if (message_len && ecc_len) {
            message = *message_len;
            ecc = *ecc_len;
            block = message + ecc;
        } else {
            throw py::value_error("Could not determine block_size and ecc_len");
        }

        if (ecc % 2 != 0)
            throw py::value_error("ecc_len must be a multiple of 2: " + std::to_string(ecc));

        return {block / sizeof(uint16_t), ecc / sizeof(uint16_t)};
    }

    template<size_t vec_size, typename T, typename U>
    inline void copy_blocks(const T *src, size_t src_size, U *dst) const {
        for (size_t j = 0; j < vec_size; ++j)
            for (size_t i = 0; i < src_size; ++i)
                dst[i + j * block_size] = src[i + j * src_size];
    }

    template<size_t vec_size, typename T, typename U>
    inline void copy_blocks_transposed(const T *src, size_t src_size, U *dst, size_t dst_block_size) const {
        for (size_t j = 0; j < vec_size; ++j)
            for (size_t i = 0; i < src_size; ++i)
                dst[i + j * dst_block_size] = src[i * vec_size + j];
    }

    template<size_t vec_size, typename T, typename U>
    inline void copy_msg_transposed(const T *src, size_t src_size, U *dst) const {
        for (size_t j = 0; j < vec_size; ++j)
            for (size_t i = 0; i < src_size; ++i)
                dst[i * vec_size + j] = src[i + j * src_size];
    }

    template<size_t vec_size, typename T>
    inline void encode_block_inline(T const& rs, const uint16_t src[], uint32_t temp[], uint16_t dst[]) {
        copy_blocks<vec_size>(&src[0], msg_size, &dst[0]);

        copy_msg_transposed<vec_size>(&src[0], msg_size, &temp[0]);
        std::fill_n(&temp[msg_size * vec_size], ecc_len * vec_size, 0);
        // std::memset(&temp[msg_size * vec_size], 0, ecc_len * vec_size * sizeof(uint32_t));

        rs.encode(&temp[0]);
        copy_blocks_transposed<vec_size>(&temp[0], ecc_len, &dst[msg_size], block_size);
    }

    template<size_t vec_size, typename T>
    inline void encode_block_external(T const& rs, const uint16_t src[], uint32_t temp[], uint16_t dst[]) {
        copy_msg_transposed<vec_size>(&src[0], msg_size, &temp[0]);
        std::fill_n(&temp[msg_size * vec_size], ecc_len * vec_size, 0);
        // std::memset(&temp[msg_size * vec_size], 0, ecc_len * vec_size * sizeof(uint32_t));

        rs.encode(&temp[0]);
        copy_blocks_transposed<vec_size>(&temp[0], ecc_len, &dst[0], ecc_len);
    }

    inline size_t encode_blocks_inline(const uint16_t src[], size_t full_blocks, uint32_t temp[], uint16_t dst[]) {
        size_t block = 0;

        if (simd_x16) {
            while (full_blocks - block >= 16) {
                encode_block_inline<16>(rs16v16, &src[block * msg_size], &temp[0], &dst[block * block_size]);
                block += 16;
            }
        }
        if (simd_x8) {
            while (full_blocks - block >= 8) {
                encode_block_inline<8>(rs16v8, &src[block * msg_size], &temp[0], &dst[block * block_size]);
                block += 8;
            }
        }
        if (simd_x4) {
            while (full_blocks - block >= 4) {
                encode_block_inline<4>(rs16v4, &src[block * msg_size], &temp[0], &dst[block * block_size]);
                block += 4;
            }
        }
        // TODO: use SIMD for all iterations
        while (block < full_blocks) {
            encode_block_inline<1>(rs16, &src[block * msg_size], &temp[0], &dst[block * block_size]);
            block += 1;
        }
        // TODO: remainder

        return block;
    }

    inline size_t encode_blocks_external(const uint16_t src[], size_t full_blocks, uint32_t temp[], uint16_t dst[]) {
        size_t block = 0;

        if (simd_x16) {
            while (full_blocks - block >= 16) {
                encode_block_external<16>(rs16v16, &src[block * msg_size], &temp[0], &dst[block * ecc_len]);
                block += 16;
            }
        }
        if (simd_x8) {
            while (full_blocks - block >= 8) {
                encode_block_external<8>(rs16v8, &src[block * msg_size], &temp[0], &dst[block * ecc_len]);
                block += 8;
            }
        }
        if (simd_x4) {
            while (full_blocks - block >= 4) {
                encode_block_external<4>(rs16v4, &src[block * msg_size], &temp[0], &dst[block * ecc_len]);
                block += 4;
            }
        }
        // TODO: use SIMD for all iterations
        while (block < full_blocks) {
            encode_block_external<1>(rs16, &src[block * msg_size], &temp[0], &dst[block * ecc_len]);
            block += 1;
        }
        // TODO: remainder

        return block;
    }
};
