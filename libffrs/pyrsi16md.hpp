/**************************************************************************
 * PyRSi16md.hpp
 *
 * Copyright 2026 Gabriel Machado
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
        size_t block_len;
        size_t ecc_len;
    };

    PyGFi16 gf;
    RSi16vImpl<uint32_t> rs16;
    RSi16v<4> rs16v4;
    RSi16v<8> rs16v8;
    RSi16v<16> rs16v16;
    bool simd_x4;
    bool simd_x8;
    bool simd_x16;

    inline PyRSi16md(
        rs_data&& args,
        uint32_t primitive,
        size_t interleave,
        bool simd_x4,
        bool simd_x8,
        bool simd_x16
    ):
        gf(primitive),
        rs16(gf, args.block_len, args.ecc_len),
        rs16v4(gf, args.block_len, args.ecc_len),
        rs16v8(gf, args.block_len, args.ecc_len),
        rs16v16(gf, args.block_len, args.ecc_len),
        simd_x4(simd_x4),
        simd_x8(simd_x8),
        simd_x16(simd_x16),
        block_len(args.block_len),
        message_len(args.block_len - args.ecc_len),
        ecc_len(args.ecc_len),
        interleave(interleave),
        chunk_len(message_len * interleave)
    { }

public:
    const size_t block_len;
    const size_t message_len;
    const size_t ecc_len;
    const size_t interleave;
    const size_t chunk_len;
    static constexpr size_t max_vec_size = 16;
    static constexpr size_t vec_align = max_vec_size * sizeof(uint32_t);

    inline PyRSi16md(
            std::optional<size_t> block_len,
            std::optional<uint16_t> message_len,
            std::optional<uint16_t> ecc_len,
            uint32_t primitive,
            size_t interleave,
            std::optional<bool> simd_x4,
            std::optional<bool> simd_x8,
            std::optional<bool> simd_x16
    ):
        PyRSi16md(
            _get_args(block_len, message_len, ecc_len),
            primitive,
            interleave,
            simd_x4.value_or(__builtin_cpu_supports("sse2")),
            simd_x8.value_or(__builtin_cpu_supports("avx2")),
            simd_x16.value_or(__builtin_cpu_supports("avx512f"))
        )
    { }

    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        py_assert(buf.size == message_len, std::to_string(buf.size));

        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_len]);
        auto output = py::bytearray(nullptr, ecc_len * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

        encode_block<1>(rs16, &buf[0], &temp[0], &output_data[0]);

        return output;
    }

    inline void repair_chunk(uint16_t message[], uint16_t ecc[], std::vector<size_t> const& error_pos) {
        auto error_pos_rbo = std::vector<size_t>(error_pos.size());
        for (size_t i = 0; i < error_pos.size(); ++i)
            error_pos_rbo[i] = rs16.rbo(error_pos[i]);

        if (simd_x16)
            _repair_chunk<16>(rs16v16, &message[0], &ecc[0], error_pos_rbo);
        else if (simd_x8)
            _repair_chunk<8>(rs16v8, &message[0], &ecc[0], error_pos_rbo);
        else if (simd_x4)
            _repair_chunk<4>(rs16v4, &message[0], &ecc[0], error_pos_rbo);
        else
            _repair_chunk<1>(rs16, &message[0], &ecc[0], error_pos_rbo);
    }

    template<size_t vec_size, typename T>
    inline void _repair_chunk(T const& rs, uint16_t message[], uint16_t ecc[], std::vector<size_t> const& error_pos_rbo) {
        // src_size = message_len * interleave
        // dst_size = ecc_len * interleave
        // block_len = message_len + ecc_len
        // chunk_size = block_len * interleave
        // temp_size = block_len * vec_size

        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_len * vec_size]);
        auto buf = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_len * vec_size]);

        size_t vec_cols = interleave / vec_size;
        for (size_t i = 0; i < vec_cols; ++i) {
            copy_stride(&message[i * vec_size], interleave, &buf[0], vec_size, vec_size, message_len);
            copy_transposed(&ecc[i * vec_size * ecc_len], ecc_len, &buf[message_len * vec_size], vec_size);

            rs.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp[0]);

            copy_stride(&buf[0], vec_size, &message[i * vec_size], interleave, vec_size, message_len);
            copy_transposed(&buf[message_len * vec_size], vec_size, &ecc[i * vec_size * ecc_len], ecc_len);
        }

        size_t encoded_cols = vec_cols * vec_size;
        if (encoded_cols < interleave) {
            size_t remaining_cols = interleave - encoded_cols;
            copy_stride(&message[encoded_cols], interleave, &buf[0], vec_size, remaining_cols, message_len);
            copy_transposed(&ecc[encoded_cols * ecc_len], ecc_len, ecc_len, &buf[message_len * vec_size], vec_size, remaining_cols);

            rs.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp[0]);

            copy_stride(&buf[0], vec_size, &message[encoded_cols], interleave, remaining_cols, message_len);
            copy_transposed(&buf[message_len * vec_size], vec_size, remaining_cols, &ecc[encoded_cols * ecc_len], ecc_len, ecc_len);
        }
    }

    inline void py_repair(buffer_rw<uint16_t> message, buffer_rw<uint16_t> ecc, std::vector<size_t> const& error_pos) {
        py_assert(message.size == message_len, std::to_string(message.size));
        py_assert(ecc.size == ecc_len, std::to_string(ecc.size));
        py_assert(error_pos.size() <= ecc_len);

        if (error_pos.empty())
            return;

        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_len]);
        auto buf = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_len]);

        std::copy_n(&message[0], message_len, &buf[0]);
        std::copy_n(&ecc[0], ecc_len, &buf[message_len]);

        auto error_pos_rbo = std::vector<size_t>(error_pos.size());
        for (size_t i = 0; i < error_pos.size(); ++i)
            error_pos_rbo[i] = rs16.rbo(error_pos[i]);

        rs16.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp[0]);

        std::copy_n(&buf[0], message_len, &message[0]);
        std::copy_n(&buf[message_len], ecc_len, &ecc[0]);
    }

    inline void encode_full_blocks(const uint16_t src[], size_t full_blocks, size_t remainder, uint16_t dst[]) {
        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_len * max_vec_size]);

        size_t block = 0;
        if (simd_x16) {
            while (full_blocks - block >= 16) {
                encode_block<16>(rs16v16, &src[block * message_len], &temp[0], &dst[block * ecc_len]);
                block += 16;
            }
        } else if (simd_x8) {
            while (full_blocks - block >= 8) {
                encode_block<8>(rs16v8, &src[block * message_len], &temp[0], &dst[block * ecc_len]);
                block += 8;
            }
        } else if (simd_x4) {
            while (full_blocks - block >= 4) {
                encode_block<4>(rs16v4, &src[block * message_len], &temp[0], &dst[block * ecc_len]);
                block += 4;
            }
        } else {
            while (block < full_blocks) {
                encode_block<1>(rs16, &src[block * message_len], &temp[0], &dst[block * ecc_len]);
                block += 1;
            }
        }

        remainder += (full_blocks - block) * message_len;

        if (remainder) {
            if (simd_x16) {
                encode_block<16>(rs16v16, &src[block * message_len], remainder, &temp[0], &dst[block * ecc_len]);
            } else if (simd_x8) {
                encode_block<8>(rs16v8, &src[block * message_len], remainder, &temp[0], &dst[block * ecc_len]);
            } else if (simd_x4) {
                encode_block<4>(rs16v4, &src[block * message_len], remainder, &temp[0], &dst[block * ecc_len]);
            } else {
                encode_block<1>(rs16, &src[full_blocks * message_len], remainder, &temp[0], &dst[full_blocks * ecc_len]);
            }
        }
    }

    inline py::bytearray py_encode_blocks(buffer_ro<uint16_t> buf) {
        if (buf.size == 0 || block_len == 0 || block_len <= ecc_len)
            return {};

        size_t full_blocks = buf.size / message_len;
        size_t output_size = full_blocks * ecc_len;

        // Last block will be smaller if input size is not divisible by message_len
        size_t input_remainder = buf.size - full_blocks * message_len;
        if (input_remainder > 0)
            output_size += ecc_len;

        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_len * max_vec_size]);
        auto output = py::bytearray(nullptr, output_size * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

        encode_full_blocks(&buf[0], full_blocks, input_remainder, &output_data[0]);

        return output;
    }

    inline void synd_blocks(const uint16_t msg[], const uint16_t ecc[], size_t count, uint32_t temp[], uint32_t synds[]) {
        // len(temp) == block_len
        // len(synds) == count * ecc_len
        for (size_t i = 0; i < count; ++i) {
            std::copy_n(&msg[i * message_len], message_len, &temp[0]);
            std::copy_n(&ecc[i * ecc_len], ecc_len, &temp[message_len]);
            rs16.ntt(&temp[0]);
            std::copy_n(&temp[0], ecc_len, &synds[i * ecc_len]);
        }
    }

    inline void synd_chunk(const uint16_t msg[], const uint16_t ecc[], uint32_t temp[], uint32_t synds[]) {
        // len(temp) == block_len
        // len(synds) == count * ecc_len

        const size_t vec_size = 1;
        for (size_t i = 0; i < interleave; ++i) {
            // std::copy_n(&msg[i * message_len], message_len, &temp[0]);
            copy_stride(&msg[i * vec_size], interleave, &temp[0], vec_size, vec_size, message_len);
            // std::copy_n(&ecc[i * ecc_len], ecc_len, &temp[message_len]);
            copy_transposed(&ecc[i * ecc_len], ecc_len, &temp[message_len * vec_size], vec_size);
            rs16.ntt(&temp[0]);
            // std::copy_n(&temp[0], ecc_len, &synds[i * ecc_len]);
            copy_transposed(&temp[0], vec_size, &synds[i * ecc_len * vec_size], ecc_len);
        }
    }

    template<typename Src, typename Dst>
    inline void encode_chunk(const Src src[], Dst dst[]) {
        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{vec_align}) uint32_t[block_len * max_vec_size]);

        if (simd_x16)
            _encode_chunk<16>(rs16v16, &src[0], &temp[0], &dst[0]);
        else if (simd_x8)
            _encode_chunk<8>(rs16v8, &src[0], &temp[0], &dst[0]);
        else if (simd_x4)
            _encode_chunk<4>(rs16v4, &src[0], &temp[0], &dst[0]);
        else
            _encode_chunk<1>(rs16, &src[0], &temp[0], &dst[0]);
    }

    inline py::bytearray py_encode_chunk(buffer_ro<uint16_t> buf) {
        py_assert(buf.size == message_len * interleave, std::to_string(buf.size));

        auto output = py::bytearray(nullptr, ecc_len * interleave * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));
        encode_chunk(&buf[0], &output_data[0]);
        return output;
    }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyRSi16md>(m, "RSi16md")
            .def_property_readonly("ecc_len", [](PyRSi16md& self) { return self.ecc_len; })
            .def_property_readonly("ecc_size", [](PyRSi16md& self) { return self.ecc_len * sizeof(uint16_t); })
            .def_property_readonly("block_len", [](PyRSi16md& self) { return self.block_len; })
            .def_property_readonly("block_size", [](PyRSi16md& self) { return self.block_len * sizeof(uint16_t); })
            .def_property_readonly("message_len", [](PyRSi16md& self) { return (self.block_len - self.ecc_len); })
            .def_property_readonly("message_size", [](PyRSi16md& self) { return self.message_len * sizeof(uint16_t); })
            .def_property_readonly("root", [](PyRSi16md& self) { return self.rs16.root; })
            .def_property_readonly("gf", [](PyRSi16md& self) -> auto const& { return self.gf; })

            .def_property("interleave",
                [](PyRSi16md& self) { return self.interleave; },
                [](PyRSi16md& self, size_t value) {
                    const_cast<size_t&>(self.interleave) = value;
                },
                R"(Number of interleaved codewords for block encoding)"
            )
            .def_property_readonly("chunk_len", [](PyRSi16md& self) { return self.chunk_len; })
            .def_property_readonly("chunk_size", [](PyRSi16md& self) { return self.chunk_len * sizeof(uint16_t); })

            .def_property("simd_x4",
                [](PyRSi16md& self) { return self.simd_x4; },
                [](PyRSi16md& self, bool value) { self.simd_x4 = value; },
                R"(Enable SIMD x4 encoding)"
            )
            .def_property("simd_x8",
                [](PyRSi16md& self) { return self.simd_x8; },
                [](PyRSi16md& self, bool value) { self.simd_x8 = value; },
                R"(Enable SIMD x8 encoding)"
            )
            .def_property("simd_x16",
                [](PyRSi16md& self) { return self.simd_x16; },
                [](PyRSi16md& self, bool value) { self.simd_x16 = value; },
                R"(Enable SIMD x16 encoding)"
            )

            .def(py::init<
                    std::optional<size_t>,    // block_size
                    std::optional<uint16_t>,  // message_size
                    std::optional<uint16_t>,  // ecc_size
                    uint32_t,  // primitive
                    size_t,    // interleave
                    std::optional<bool>,
                    std::optional<bool>,
                    std::optional<bool>
                >(),
                R"(Instantiate a Reed-Solomon encoder with the given configuration)",
                "block_len"_a = py::none(),
                "message_len"_a = py::none(),
                "ecc_len"_a = py::none(),
                "primitive"_a = 3,
                "interleave"_a = 1,
                "simd_x4"_a = py::none(),
                "simd_x8"_a = py::none(),
                "simd_x16"_a = py::none()
            )

            .def("__sizeof__", [](PyRSi16md& self) { return sizeof(self); })

            .def("encode", cast_args(&PyRSi16md::py_encode),
                R"(Systematic encode)",
                "buffer"_a)

            .def("encode_blocks", cast_args(&PyRSi16md::py_encode_blocks), R"(
                Encode the input buffer in blocks, storing the parity bytes right next to their corresponding block
                )",
                "buffer"_a)

            .def("encode_chunk", cast_args(&PyRSi16md::py_encode_chunk),
                R"(Encode a chunk of interleaved codewords)",
                "buffer"_a)

            .def("repair", cast_args(&PyRSi16md::py_repair),
                R"(Repair a block with the given error locations)",
                "message"_a,
                "ecc"_a,
                "error_pos"_a)

            // .def("decode", cast_args(&PyRSi16md::py_decode),
            //     R"(Systematic decode)",
            //     "buffer"_a)

            .doc() = R"(Reed-Solomon coding over :math:`GF(2^8)`)";
    }

private:
    inline rs_data _get_args(
            std::optional<size_t> block_len,
            std::optional<uint16_t> message_len,
            std::optional<uint16_t> ecc_len) {

        size_t block, message, ecc;
        if (ecc_len && block_len) {
            ecc = *ecc_len;
            block = *block_len;
            message = block - ecc;

            if (message_len and *message_len != block - ecc)
                throw py::value_error("block_len = message_len + ecc_len");
        } else if (message_len && block_len) {
            if (*message_len >= *block_len)
                throw py::value_error("block_len must be greater than message_len");

            block = *block_len;
            message = *message_len;
            ecc = block - message;
        } else if (message_len && ecc_len) {
            message = *message_len;
            ecc = *ecc_len;
            block = message + ecc;
        } else {
            throw py::value_error("Could not determine block_len and ecc_len");
        }

        if (ecc % 2 != 0)
            throw py::value_error("ecc_len must be a multiple of 2: " + std::to_string(ecc));

        return {block, ecc};
    }

    template<typename Src, typename Dst>
    inline void copy_transposed(const Src src[], size_t src_cols, Dst dst[], size_t dst_cols) const {
        for (size_t j = 0; j < src_cols; ++j)
            for (size_t i = 0; i < dst_cols; ++i)
                dst[i + j * dst_cols] = src[i * src_cols + j];
    }

    template<typename Src, typename Dst>
    inline void copy_transposed(const Src src[], size_t src_stride, size_t src_cols, Dst dst[], size_t dst_stride, size_t dst_cols) const {
        for (size_t j = 0; j < src_cols; ++j)
            for (size_t i = 0; i < dst_cols; ++i)
                dst[i + j * dst_stride] = src[i * src_stride + j];
    }

    template<typename Src, typename Dst>
    inline void copy_stride(const Src *src, size_t src_stride, Dst *dst, size_t dst_stride, size_t cols, size_t count) const {
        size_t src_i = 0;
        size_t dst_i = 0;
        for (size_t i = 0; i < count; ++i) {
            std::copy_n(&src[src_i], cols, &dst[dst_i]);
            // std::memcpy(&dst[dst_i], &src[src_i], cols * sizeof(T));
            src_i += src_stride;
            dst_i += dst_stride;
        }

        // for (size_t i = 0; i < count; ++i) {
        //     // std::copy_n(src + i * src_stride, n, dst + i * dst_stride);
        //     std::memcpy(dst + i * dst_stride, src + i * src_stride, cols * sizeof(T));
        // }
    }

    template<size_t vec_size, typename RS, typename Src, typename Dst>
    inline void _encode_chunk(RS const& rs, const Src src[], uint32_t temp[], Dst dst[]) {
        // src_size = message_len * interleave
        // dst_size = ecc_len * interleave
        // block_len = message_len + ecc_len
        // chunk_size = block_len * interleave
        // temp_size = block_len * vec_size
        size_t vec_cols = interleave / vec_size;
        for (size_t i = 0; i < vec_cols; ++i) {
            copy_stride(&src[i * vec_size], interleave, &temp[0], vec_size, vec_size, message_len);
            std::fill_n(&temp[message_len * vec_size], ecc_len * vec_size, 0);
            // std::memset(&temp[message_len * vec_size], 0, ecc_len * vec_size * sizeof(uint32_t));

            rs.encode(&temp[0]);
            // std::copy_n(&temp[0], ecc_len * vec_size, &dst[i * vec_size * ecc_len]);
            copy_transposed(&temp[0], vec_size, &dst[i * ecc_len * vec_size], ecc_len);
            // copy_stride(&temp[0], vec_size, &dst[i * vec_size], interleave, vec_size, ecc_len);
        }
        size_t encoded_cols = vec_cols * vec_size;
        if (encoded_cols < interleave) {
            size_t remaining_cols = interleave - encoded_cols;
            copy_stride(&src[encoded_cols], interleave, &temp[0], vec_size, remaining_cols, message_len);
            std::fill_n(&temp[message_len * vec_size], ecc_len * vec_size, 0);
            // std::memset(&temp[message_len * vec_size], 0, ecc_len * vec_size * sizeof(uint32_t));

            rs.encode(&temp[0]);
            // copy_stride(&temp[0], vec_size, &dst[encoded_cols], remaining_cols, remaining_cols, ecc_len);
            copy_transposed(&temp[0], vec_size, remaining_cols, &dst[encoded_cols * ecc_len], ecc_len, ecc_len);
            // copy_stride(&temp[0], vec_size, &dst[encoded_cols], interleave, remaining_cols, ecc_len);
        }
    }

    template<size_t vec_size, typename T>
    inline void encode_block(T const& rs, const uint16_t src[], uint32_t temp[], uint16_t dst[]) {
        copy_transposed(&src[0], message_len, &temp[0], vec_size);
        std::fill_n(&temp[message_len * vec_size], ecc_len * vec_size, 0);
        // std::memset(&temp[message_len * vec_size], 0, ecc_len * vec_size * sizeof(uint32_t));

        rs.encode(&temp[0]);
        copy_transposed(&temp[0], vec_size, &dst[0], ecc_len);
    }

    template<size_t vec_size, typename T>
    inline void encode_block(T const& rs, const uint16_t src[], size_t src_size, uint32_t temp[], uint16_t dst[]) {
        std::fill_n(&temp[0], block_len * vec_size, 0);

        auto cols = src_size / message_len;
        copy_transposed(&src[0], message_len, message_len, &temp[0], vec_size, cols);

        size_t remainder = src_size - cols * message_len;
        if (remainder) {
            copy_transposed(&src[src_size - remainder], remainder, remainder, &temp[cols], vec_size, 1);
            ++cols;
        }

        rs.encode(&temp[0]);
        copy_transposed(&temp[0], vec_size, cols, &dst[0], ecc_len, ecc_len);
    }

    template<typename T>
    inline void _print_table(T const& table, size_t rows, size_t cols) const {
        py::print("Table:", rows, "x", cols);
        for (size_t i = 0; i < rows; ++i) {
            std::string row;
            for (size_t j = 0; j < cols; ++j) {
                auto val = table[i * cols + j];
                if (val < 10)
                    row += " ";
                if (val < 100)
                    row += " ";
                row += std::to_string(val) + " ";
            }
            py::print(row);
        }
    }
};
