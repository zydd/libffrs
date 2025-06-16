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

    PyGFi16 gf;
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

    inline PyRSi16md(
        rs_data&& args,
        uint32_t primitive,
        std::optional<bool> simd_x16,
        std::optional<bool> simd_x8,
        std::optional<bool> simd_x4
    ):
        gf(primitive),
        rs16(gf, args.block_size, args.ecc_len),
        rs16v4(gf, args.block_size, args.ecc_len),
        rs16v8(gf, args.block_size, args.ecc_len),
        rs16v16(gf, args.block_size, args.ecc_len),
        block_size(args.block_size),
        msg_size(args.block_size - args.ecc_len),
        ecc_len(args.ecc_len),
        simd_x16(simd_x16.value_or(__builtin_cpu_supports("avx512f"))),
        simd_x8(simd_x8.value_or(__builtin_cpu_supports("avx2"))),
        simd_x4(simd_x4.value_or(__builtin_cpu_supports("sse2")))
    { }

public:
    static constexpr size_t max_vec_size = 16;
    static constexpr size_t vec_align = max_vec_size * sizeof(uint32_t);

    inline PyRSi16md(
            std::optional<size_t> block_size,
            std::optional<uint16_t> message_len,
            std::optional<uint16_t> ecc_len,
            uint32_t primitive,
            std::optional<bool> simd_x16,
            std::optional<bool> simd_x8,
            std::optional<bool> simd_x4
    ):
        PyRSi16md(
            _get_args(block_size, message_len, ecc_len),
            primitive,
            simd_x16,
            simd_x8,
            simd_x4
        )
    { }

    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        constexpr size_t vec_size = 1;
        py_assert(buf.size == msg_size * vec_size, std::to_string(buf.size));

        size_t output_size = ecc_len;

        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_size * vec_size]);
        auto output = py::bytearray(nullptr, output_size * vec_size * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

        encode_block<vec_size>(rs16, &buf[0], &temp[0], &output_data[0]);

        return output;
    }

    inline py::bytearray py_encode_blocks(buffer_ro<uint16_t> buf) {
        if (buf.size == 0 || block_size == 0 || block_size <= ecc_len)
            return {};

        size_t full_blocks = buf.size / msg_size;
        size_t output_size = full_blocks * ecc_len;

        // Last block will be smaller if input size is not divisible by msg_size
        size_t input_remainder = buf.size - full_blocks * msg_size;
        if (input_remainder > 0)
            output_size += ecc_len;

        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_size * max_vec_size]);
        auto output = py::bytearray(nullptr, output_size * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

        size_t block = 0;
        if (simd_x16) {
            while (full_blocks - block >= 16) {
                encode_block<16>(rs16v16, &buf[block * msg_size], &temp[0], &output_data[block * ecc_len]);
                block += 16;
            }
            input_remainder = buf.size - block * msg_size;
            if (input_remainder)
                encode_block<16>(rs16v16, &buf[block * msg_size], input_remainder, &temp[0], &output_data[block * ecc_len]);
        } else if (simd_x8) {
            while (full_blocks - block >= 8) {
                encode_block<8>(rs16v8, &buf[block * msg_size], &temp[0], &output_data[block * ecc_len]);
                block += 8;
            }
            input_remainder = buf.size - block * msg_size;
            if (input_remainder)
                encode_block<8>(rs16v8, &buf[block * msg_size], input_remainder, &temp[0], &output_data[block * ecc_len]);
        } else if (simd_x4) {
            while (full_blocks - block >= 4) {
                encode_block<4>(rs16v4, &buf[block * msg_size], &temp[0], &output_data[block * ecc_len]);
                block += 4;
            }
            input_remainder = buf.size - block * msg_size;
            if (input_remainder)
                encode_block<4>(rs16v4, &buf[block * msg_size], input_remainder, &temp[0], &output_data[block * ecc_len]);
        } else {
            while (block < full_blocks) {
                encode_block<1>(rs16, &buf[block * msg_size], &temp[0], &output_data[block * ecc_len]);
                block += 1;
            }
            if (input_remainder)
                encode_block<1>(rs16, &buf[block * msg_size], input_remainder, &temp[0], &output_data[block * ecc_len]);
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

            .def_property_readonly("root",
                [](PyRSi16md& self) { return self.rs16.root; })

            .def_property_readonly("gf", [](PyRSi16md& self) -> auto const& { return self.gf; })

            .def(py::init<
                    std::optional<uint16_t>,
                    std::optional<uint16_t>,
                    std::optional<uint16_t>,
                    uint32_t,
                    std::optional<bool>,
                    std::optional<bool>,
                    std::optional<bool>>(),
                R"(Instantiate a Reed-Solomon encoder with the given configuration)",
                "block_size"_a = py::none(),
                "message_len"_a = py::none(),
                "ecc_len"_a = py::none(),
                "primitive"_a = 3,
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
    inline void copy_ecc_transposed(const T *src, U *dst, size_t dst_stride) const {
        for (size_t j = 0; j < vec_size; ++j)
            for (size_t i = 0; i < ecc_len; ++i)
                dst[i + j * dst_stride] = src[i * vec_size + j];
    }

    template<size_t vec_size, typename T, typename U>
    inline void copy_ecc_transposed(const T *src, size_t cols, U *dst, size_t dst_stride) const {
        for (size_t j = 0; j < cols; ++j)
            for (size_t i = 0; i < ecc_len; ++i)
                dst[i + j * dst_stride] = src[i * vec_size + j];
    }

    template<size_t vec_size, typename T, typename U>
    inline void copy_msg_transposed(const T *src, U *dst) const {
        for (size_t j = 0; j < vec_size; ++j)
            for (size_t i = 0; i < msg_size; ++i)
                dst[i * vec_size + j] = src[i + j * msg_size];
    }

    template<size_t vec_size, typename T, typename U>
    inline size_t copy_msg_transposed(const T *src, size_t src_size, U *dst) const {
        size_t i, j;
        for (j = 0; j < vec_size; ++j) {
            for (i = 0; i < msg_size; ++i) {
                auto pos = i + j * msg_size;
                if (pos >= src_size) [[unlikely]] {
                    if (i > 0) {
                        // Empty remaining column
                        for (; i < msg_size; ++i)
                            dst[i * vec_size + j] = 0;
                        return j + 1;
                    }
                    return j;
                }
                dst[i * vec_size + j] = src[pos];
            }
        }

        return j;
    }

    template<size_t vec_size, typename T>
    inline void encode_block(T const& rs, const uint16_t src[], uint32_t temp[], uint16_t dst[]) {
        copy_msg_transposed<vec_size>(&src[0], &temp[0]);
        std::fill_n(&temp[msg_size * vec_size], ecc_len * vec_size, 0);
        // std::memset(&temp[msg_size * vec_size], 0, ecc_len * vec_size * sizeof(uint32_t));

        rs.encode(&temp[0]);
        copy_ecc_transposed<vec_size>(&temp[0], &dst[0], ecc_len);
    }

    template<size_t vec_size, typename T>
    inline void encode_block(T const& rs, const uint16_t src[], size_t src_size, uint32_t temp[], uint16_t dst[]) {
        auto cols = copy_msg_transposed<vec_size>(&src[0], src_size, &temp[0]);
        std::fill_n(&temp[msg_size * vec_size], ecc_len * vec_size, 0);
        // std::memset(&temp[msg_size * vec_size], 0, ecc_len * vec_size * sizeof(uint32_t));

        rs.encode(&temp[0]);
        copy_ecc_transposed<vec_size>(&temp[0], cols, &dst[0], ecc_len);
    }
};
