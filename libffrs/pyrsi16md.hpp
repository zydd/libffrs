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
// #include "rsi16md_impl.hpp"
#include "rsi16md.h"

namespace py = pybind11;


class PyRSi16md {
private:
    struct rs_data {
        size_t block_len;
        size_t ecc_len;
    };

    PyGFi16 gf;
    RSi16v<1> rs16;
    RSi16v<4> rs16x4;
    RSi16v<8> rs16x8;
    RSi16v<16> rs16x16;
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
        rs16x4(gf, args.block_len, args.ecc_len),
        rs16x8(gf, args.block_len, args.ecc_len),
        rs16x16(gf, args.block_len, args.ecc_len),
        simd_x4(simd_x4),
        simd_x8(simd_x8),
        simd_x16(simd_x16),
        block_len(args.block_len),
        message_len(args.block_len - args.ecc_len),
        ecc_len(args.ecc_len),
        interleave(interleave),
        interleaved_block_len(block_len * interleave),
        interleaved_message_len(message_len * interleave),
        interleaved_ecc_len(ecc_len * interleave)
    {
        if (simd_x16)
            vec_align = 16 * sizeof(GFT);
        else if (simd_x8)
            vec_align = 8 * sizeof(GFT);
        else if (simd_x4)
            vec_align = 4 * sizeof(GFT);
        else
            vec_align = sizeof(GFT);
    }

public:
    const size_t block_len;
    const size_t message_len;
    const size_t ecc_len;
    const size_t interleave;
    const size_t interleaved_block_len;
    const size_t interleaved_message_len;
    const size_t interleaved_ecc_len;
    size_t vec_align;

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
            simd_x16.value_or(__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"))
        )
    { }

    template<typename Src, typename Dst>
    inline void encode_blocks(const Src src[], size_t full_blocks, size_t remainder, Dst dst[]) {
        size_t block = 0;
        _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            auto temp = new_aligned<GFT>(block_len * SIMD_W, SIMD_W * sizeof(GFT));
            while (full_blocks - block >= SIMD_W) {
                encode_block<SIMD_W>(rs, &src[block * message_len], &temp[0], &dst[block * ecc_len]);
                block += SIMD_W;
            }

            remainder += (full_blocks - block) * message_len;
            if (remainder) {
                encode_block<SIMD_W>(rs, &src[block * message_len], remainder, &temp[0], &dst[block * ecc_len]);
            }
        });
    }

    template<typename Src, typename Dst>
    inline void encode_interleaved(const Src src[], Dst dst[]) {
        _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            auto temp = new_aligned<GFT>(block_len * SIMD_W, SIMD_W * sizeof(GFT));
            _encode_interleaved<SIMD_W>(rs, &src[0], &temp[0], &dst[0]);
        });
    }

    template<typename Msg, typename Ecc>
    inline void repair_interleaved(Msg message[], Ecc ecc[], size_t col_start, size_t col_count) {
        _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            _repair_interleaved<SIMD_W>(rs, &message[0], &ecc[0], col_start, col_count);
        });
    }

    template<typename Msg, typename Ecc>
    inline void repair_interleaved(Msg message[], Ecc ecc[], size_t col_start, size_t col_count, std::vector<size_t> const& error_pos) {
        auto error_pos_rbo = std::vector<size_t>(error_pos.size());
        for (size_t i = 0; i < error_pos.size(); ++i)
            error_pos_rbo[i] = rs16.rbo(error_pos[i]);

        _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            _repair_interleaved<SIMD_W>(rs, &message[0], &ecc[0], col_start, col_count, error_pos_rbo);
        });
    }

    inline void repair_block(uint32_t block[], std::vector<size_t> const& error_pos, uint32_t temp2[]) const {
        auto error_pos_rbo = std::vector<size_t>(error_pos.size());
        for (size_t i = 0; i < error_pos.size(); ++i)
            error_pos_rbo[i] = rs16.rbo(error_pos[i]);

        py_assert(block_len == rs16.ntt_len, "block len != NTT len: " + std::to_string(block_len) + " != " + std::to_string(rs16.ntt_len));
        rs16.repair(&block[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp2[0]);
    }

    template<typename Src>
    inline void synd_blocks(const Src msg[], const uint16_t ecc[], size_t count, uint32_t temp[], uint32_t synds[]) {
        // len(temp) == block_len
        // len(synds) == count * ecc_len
        for (size_t i = 0; i < count; ++i) {
            std::copy_n(&msg[i * message_len], message_len, &temp[0]);
            std::copy_n(&ecc[i * ecc_len], ecc_len, &temp[message_len]);
            rs16.pntt(&temp[0]);
            std::copy_n(&temp[0], ecc_len, &synds[i * ecc_len]);
        }
    }

    inline void synd_interleaved(const uint16_t msg[], const uint16_t ecc[], uint32_t temp[], uint32_t synds[]) {
        // len(temp) == block_len
        // len(synds) == count * ecc_len

        const size_t SIMD_W = 1;
        for (size_t i = 0; i < interleave; ++i) {
            // std::copy_n(&msg[i * message_len], message_len, &temp[0]);
            copy_stride(&msg[i * SIMD_W], interleave, &temp[0], SIMD_W, SIMD_W, message_len);
            // std::copy_n(&ecc[i * ecc_len], ecc_len, &temp[message_len]);
            copy_transposed(&ecc[i * ecc_len], ecc_len, &temp[message_len * SIMD_W], SIMD_W);
            rs16.pntt(&temp[0]);
            // std::copy_n(&temp[0], ecc_len, &synds[i * ecc_len]);
            copy_transposed(&temp[0], SIMD_W, &synds[i * ecc_len * SIMD_W], ecc_len);
        }
    }

    inline void ecc_mix(uint32_t ecc[]) {
        rs16.ecc_mix(&ecc[0]);
    }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyRSi16md>(m, "RSi16md")
            .def_property_readonly("ecc_len", [](PyRSi16md& self) { return self.interleaved_ecc_len; })
            .def_property_readonly("ecc_size", [](PyRSi16md& self) { return self.interleaved_ecc_len * sizeof(uint16_t); })
            .def_property_readonly("block_len", [](PyRSi16md& self) { return self.interleaved_block_len; })
            .def_property_readonly("block_size", [](PyRSi16md& self) { return self.interleaved_block_len * sizeof(uint16_t); })
            .def_property_readonly("message_len", [](PyRSi16md& self) { return self.interleaved_message_len; })
            .def_property_readonly("message_size", [](PyRSi16md& self) { return self.interleaved_message_len * sizeof(uint16_t); })
            .def_property_readonly("ntt_len", [](PyRSi16md& self) { return self.rs16.ntt_len; })
            .def_property_readonly("ntt_size", [](PyRSi16md& self) { return self.rs16.ntt_len * sizeof(uint16_t); })
            .def_property_readonly("root", [](PyRSi16md& self) { return self.rs16.root; })
            .def_property_readonly("gf", [](PyRSi16md& self) -> auto const& { return self.gf; })
            .def_property_readonly("interleave", [](PyRSi16md& self) { return self.interleave; },
                R"(Number of interleaved codewords for block encoding)")

            .def_property_readonly("simd_x4",
                [](PyRSi16md& self) { return self.simd_x4; },
                R"(Enable SIMD x4 encoding)"
            )
            .def_property_readonly("simd_x8",
                [](PyRSi16md& self) { return self.simd_x8; },
                R"(Enable SIMD x8 encoding)"
            )
            .def_property_readonly("simd_x16",
                [](PyRSi16md& self) { return self.simd_x16; },
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

            .def("synd", cast_args(&PyRSi16md::py_synd),
                R"(Calculate syndromes for the given message and ecc buffers)",
                "message"_a,
                "ecc"_a)

            .def("repair", cast_args(&PyRSi16md::py_repair),
                R"(Repair a block with the given error locations)",
                "message"_a,
                "ecc"_a,
                "error_pos"_a = py::none())

            .doc() = R"(Reed-Solomon coding over :math:`GF(65537)`)";
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

    template<size_t SIMD_W, typename RS, typename Src, typename Dst>
    inline void _encode_interleaved(RS const& rs, const Src src[], uint32_t temp[], Dst dst[]) {
        // src_size = message_len * interleave
        // dst_size = ecc_len * interleave
        // block_len = message_len + ecc_len
        // chunk_size = block_len * interleave
        // temp_size = block_len * SIMD_W
        size_t vec_cols = interleave / SIMD_W;
        for (size_t i = 0; i < vec_cols; ++i) {
            copy_stride(&src[i * SIMD_W], interleave, &temp[0], SIMD_W, SIMD_W, message_len);
            std::fill_n(&temp[message_len * SIMD_W], ecc_len * SIMD_W, 0);
            // std::memset(&temp[message_len * SIMD_W], 0, ecc_len * SIMD_W * sizeof(GFT));

            rs.encode(&temp[0]);
            // std::copy_n(&temp[0], ecc_len * SIMD_W, &dst[i * SIMD_W * ecc_len]);

            // Sequential ecc
            // copy_transposed(&temp[0], SIMD_W, &dst[i * ecc_len * SIMD_W], ecc_len);

            // Interleaved ecc
            copy_stride(&temp[0], SIMD_W, &dst[i * SIMD_W], interleave, SIMD_W, ecc_len);
        }
        size_t encoded_cols = vec_cols * SIMD_W;
        if (encoded_cols < interleave) {
            size_t remaining_cols = interleave - encoded_cols;
            copy_stride(&src[encoded_cols], interleave, &temp[0], SIMD_W, remaining_cols, message_len);
            std::fill_n(&temp[message_len * SIMD_W], ecc_len * SIMD_W, 0);
            // std::memset(&temp[message_len * SIMD_W], 0, ecc_len * SIMD_W * sizeof(GFT));

            rs.encode(&temp[0]);
            // copy_stride(&temp[0], SIMD_W, &dst[encoded_cols], remaining_cols, remaining_cols, ecc_len);

            // Sequential ecc
            // copy_transposed(&temp[0], SIMD_W, remaining_cols, &dst[encoded_cols * ecc_len], ecc_len, ecc_len);

            // Interleaved ecc
            copy_stride(&temp[0], SIMD_W, &dst[encoded_cols], interleave, remaining_cols, ecc_len);
        }
    }

    template<size_t SIMD_W, typename T, typename Src, typename Dst>
    inline void encode_block(T const& rs, const Src src[], uint32_t temp[], Dst dst[]) {
        copy_transposed(&src[0], message_len, &temp[0], SIMD_W);
        std::fill_n(&temp[message_len * SIMD_W], ecc_len * SIMD_W, 0);
        // std::memset(&temp[message_len * SIMD_W], 0, ecc_len * SIMD_W * sizeof(GFT));

        rs.encode(&temp[0]);
        copy_transposed(&temp[0], SIMD_W, &dst[0], ecc_len);
    }

    template<size_t SIMD_W, typename T, typename Src, typename Dst>
    inline void encode_block(T const& rs, const Src src[], size_t src_size, uint32_t temp[], Dst dst[]) {
        std::fill_n(&temp[0], block_len * SIMD_W, 0);

        auto cols = src_size / message_len;
        copy_transposed(&src[0], message_len, message_len, &temp[0], SIMD_W, cols);

        size_t remainder = src_size - cols * message_len;
        if (remainder) {
            copy_transposed(&src[src_size - remainder], remainder, remainder, &temp[cols], SIMD_W, 1);
            ++cols;
        }

        rs.encode(&temp[0]);
        copy_transposed(&temp[0], SIMD_W, cols, &dst[0], ecc_len, ecc_len);
    }

    template<size_t SIMD_W, typename Rs, typename Msg, typename Ecc>
    inline void _repair_interleaved(Rs const& rs, Msg message[], Ecc ecc[], size_t col_start, size_t col_count) {
        // src_size = message_len * interleave
        // dst_size = ecc_len * interleave
        // block_len = message_len + ecc_len
        // chunk_size = block_len * interleave

        auto temp1_ecc6 = new_aligned<GFT>(std::max(block_len, ecc_len * 6) * SIMD_W, SIMD_W * sizeof(GFT));
        auto buf = new_aligned<GFT>(block_len * SIMD_W, SIMD_W * sizeof(GFT));

        message += col_start;

        // Sequential ecc
        // ecc += col_start * ecc_len;

        // Interleaved ecc
        ecc += col_start;

        size_t vec_cols = col_count / SIMD_W;
        for (size_t i = 0; i < vec_cols; ++i) {
            copy_stride(&message[i * SIMD_W], interleave, &buf[0], SIMD_W, SIMD_W, message_len);

            // Sequential ecc
            // copy_transposed(&ecc[i * SIMD_W * ecc_len], ecc_len, &buf[message_len * SIMD_W], SIMD_W);

            // Interleaved ecc
            copy_stride(&ecc[i * SIMD_W], interleave, &buf[message_len * SIMD_W], SIMD_W, SIMD_W, ecc_len);

            rs.repair(&buf[0], &temp1_ecc6[0]);

            copy_stride(&buf[0], SIMD_W, &message[i * SIMD_W], interleave, SIMD_W, message_len);

            // Sequential ecc
            // copy_transposed(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W * ecc_len], ecc_len);

            // Interleaved ecc
            copy_stride(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W], interleave, SIMD_W, ecc_len);
        }

        size_t encoded_cols = vec_cols * SIMD_W;
        if (encoded_cols < col_count) {
            std::fill_n(&buf[0], block_len * SIMD_W, GFT{0});

            size_t remaining_cols = col_count - encoded_cols;
            copy_stride(&message[encoded_cols], interleave, &buf[0], SIMD_W, remaining_cols, message_len);

            // Sequential ecc
            // copy_transposed(&ecc[encoded_cols * ecc_len], ecc_len, ecc_len, &buf[message_len * SIMD_W], SIMD_W, remaining_cols);

            // Interleaved ecc
            copy_stride(&ecc[encoded_cols], interleave, &buf[message_len * SIMD_W], SIMD_W, remaining_cols, ecc_len);

            rs.repair(&buf[0], &temp1_ecc6[0]);

            copy_stride(&buf[0], SIMD_W, &message[encoded_cols], interleave, remaining_cols, message_len);

            // Sequential ecc
            // copy_transposed(&buf[message_len * SIMD_W], SIMD_W, remaining_cols, &ecc[encoded_cols * ecc_len], ecc_len, ecc_len);

            // Interleaved ecc
            copy_stride(&buf[message_len * SIMD_W], SIMD_W, &ecc[encoded_cols], interleave, remaining_cols, ecc_len);
        }
    }

    template<size_t SIMD_W, typename Rs, typename Msg, typename Ecc>
    inline void _repair_interleaved(Rs const& rs, Msg message[], Ecc ecc[], size_t col_start, size_t col_count, std::vector<size_t> const& error_pos_rbo) {
        // src_size = message_len * interleave
        // dst_size = ecc_len * interleave
        // block_len = message_len + ecc_len
        // interleaved_size = block_len * interleave
        // temp2_size = block_len * SIMD_W * 2

        auto temp2 = new_aligned<GFT>(block_len * SIMD_W * 2, SIMD_W * sizeof(GFT));

        py_assert(block_len == rs.ntt_len, "block len != NTT len: " + std::to_string(block_len) + " != " + std::to_string(rs.ntt_len));
        auto buf = new_aligned<GFT>(block_len * SIMD_W, SIMD_W * sizeof(GFT));

        message += col_start;

        // Sequential ecc
        // ecc += col_start * ecc_len;

        // Interleaved ecc
        ecc += col_start;

        size_t vec_cols = col_count / SIMD_W;
        for (size_t i = 0; i < vec_cols; ++i) {
            copy_stride(&message[i * SIMD_W], interleave, &buf[0], SIMD_W, SIMD_W, message_len);

            // Sequential ecc
            // copy_transposed(&ecc[i * SIMD_W * ecc_len], ecc_len, &buf[message_len * SIMD_W], SIMD_W);

            // Interleaved ecc
            copy_stride(&ecc[i * SIMD_W], interleave, &buf[message_len * SIMD_W], SIMD_W, SIMD_W, ecc_len);

            rs.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp2[0]);

            copy_stride(&buf[0], SIMD_W, &message[i * SIMD_W], interleave, SIMD_W, message_len);

            // Sequential ecc
            // copy_transposed(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W * ecc_len], ecc_len);

            // Interleaved ecc
            copy_stride(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W], interleave, SIMD_W, ecc_len);
        }

        size_t encoded_cols = vec_cols * SIMD_W;
        if (encoded_cols < col_count) {
            size_t remaining_cols = col_count - encoded_cols;
            copy_stride(&message[encoded_cols], interleave, &buf[0], SIMD_W, remaining_cols, message_len);

            // Sequential ecc
            // copy_transposed(&ecc[encoded_cols * ecc_len], ecc_len, ecc_len, &buf[message_len * SIMD_W], SIMD_W, remaining_cols);

            // Interleaved ecc
            copy_stride(&ecc[encoded_cols], interleave, &buf[message_len * SIMD_W], SIMD_W, remaining_cols, ecc_len);

            rs.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp2[0]);

            copy_stride(&buf[0], SIMD_W, &message[encoded_cols], interleave, remaining_cols, message_len);

            // Sequential ecc
            // copy_transposed(&buf[message_len * SIMD_W], SIMD_W, remaining_cols, &ecc[encoded_cols * ecc_len], ecc_len, ecc_len);

            // Interleaved ecc
            copy_stride(&buf[message_len * SIMD_W], SIMD_W, &ecc[encoded_cols], interleave, remaining_cols, ecc_len);
        }
    }

    template<typename F>
    inline void _simd_dispatch(F&& f) {
        if (simd_x16)
            f(std::integral_constant<size_t, 16>{}, rs16x16);
        else if (simd_x8)
            f(std::integral_constant<size_t, 8>{}, rs16x8);
        else if (simd_x4)
            f(std::integral_constant<size_t, 4>{}, rs16x4);
        else
            f(std::integral_constant<size_t, 1>{}, rs16);
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

    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        if (buf.size == 0)
            return {};

        if (interleave == 1) {
            size_t full_blocks = buf.size / message_len;
            size_t output_size = full_blocks * ecc_len;

            // Last block will be smaller if input size is not divisible by message_len
            size_t input_remainder = buf.size - full_blocks * message_len;
            if (input_remainder > 0)
                output_size += ecc_len;

            auto output = py::bytearray(nullptr, output_size * sizeof(uint16_t));
            auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

            encode_blocks(&buf[0], full_blocks, input_remainder, &output_data[0]);
            return output;
        } else {
            // TODO: encode multiple interleaved chunks
            py_assert(buf.size == message_len * interleave, std::to_string(buf.size));

            auto output = py::bytearray(nullptr, ecc_len * interleave * sizeof(uint16_t));
            auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));
            encode_interleaved(&buf[0], &output_data[0]);
            return output;
        }
    }

    inline void py_repair(buffer_rw<uint16_t> message, buffer_rw<uint16_t> ecc, std::optional<std::vector<size_t>> const& error_pos) {
        if (error_pos && error_pos->empty())
            return;

        if (interleave == 1) {
            py_assert(message.size == message_len, std::to_string(message.size));
            py_assert(ecc.size == ecc_len, std::to_string(ecc.size));

            auto buf = new_aligned<GFT>(block_len, vec_align);

            std::copy_n(&message[0], message_len, &buf[0]);
            std::copy_n(&ecc[0], ecc_len, &buf[message_len]);

            if (error_pos) {
                py_assert(error_pos->size() <= ecc_len);

                auto error_pos_rbo = std::vector<size_t>(error_pos->size());
                for (size_t i = 0; i < error_pos->size(); ++i)
                    error_pos_rbo[i] = rs16.rbo((*error_pos)[i]);

                auto temp2 = new_aligned<GFT>(block_len * 2, sizeof(GFT));
                py_assert(block_len == rs16.ntt_len, "block len != NTT len: " + std::to_string(block_len) + " != " + std::to_string(rs16.ntt_len));
                rs16.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp2[0]);

            } else {
                auto temp1_ecc6 = new_aligned<GFT>(std::max(block_len, ecc_len * 6), sizeof(GFT));
                rs16.repair(&buf[0], &temp1_ecc6[0]);
            }

            std::copy_n(&buf[0], message_len, &message[0]);
            std::copy_n(&buf[message_len], ecc_len, &ecc[0]);
        } else {
            py_assert(message.size == interleaved_message_len, std::to_string(message.size));
            py_assert(ecc.size == interleaved_ecc_len, std::to_string(ecc.size));

            if (error_pos) {
                py_assert(error_pos->size() <= ecc_len);
                repair_interleaved(&message[0], &ecc[0], 0, interleave, *error_pos);
            } else {
                repair_interleaved(&message[0], &ecc[0], 0, interleave);
            }

        }
    }

    inline py::bytearray py_synd(buffer_rw<uint16_t> message, buffer_rw<uint16_t> ecc) {
        py_assert(message.size == message_len, std::to_string(message.size));
        py_assert(ecc.size == ecc_len, std::to_string(ecc.size));

        auto temp = new_aligned<GFT>(block_len, sizeof(GFT));
        auto synds = new_aligned<GFT>(ecc_len, sizeof(GFT));

        std::copy_n(&message[0], message_len, &temp[0]);
        std::copy_n(&ecc[0], ecc_len, &temp[message_len]);

        rs16.pntt(&temp[0]);
        std::copy_n(&temp[0], ecc_len, &synds[0]);

        auto output = py::bytearray(nullptr, ecc_len * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));
        std::copy_n(&synds[0], ecc_len, &output_data[0]);
        return output;
    }
};
