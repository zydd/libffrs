/**************************************************************************
 * PyRSi16.hpp
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
// #include "rsi16v_impl.hpp"
#include "rsi16v.h"

namespace py = pybind11;


class PyRSi16 {
private:
    struct rs_data {
        size_t block_len;
        size_t ecc_len;
    };

    inline PyRSi16(
        rs_data&& args,
        size_t interleave,
        uint32_t primitive,
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
        interleaved_ecc_len(ecc_len * interleave),
        repair_temp_len(rs16.ntt_len + ecc_len * 6)
    {
        if (simd_x16)
            vec_align = 16 * sizeof(GFT);
        else if (simd_x8)
            vec_align = 8 * sizeof(GFT);
        else if (simd_x4)
            vec_align = 4 * sizeof(GFT);
        else
            vec_align = sizeof(GFT);

        py_assert(
            ecc_len >= 2 && (ecc_len & (ecc_len - 1)) == 0,
            "ecc_len must be a power of 2: " + std::to_string(ecc_len)
        );
        py_assert(
            block_len % ecc_len == 0,
            "block_len must be a multiple of ecc_len. block_len: "
            + std::to_string(block_len) + " ecc_len: " + std::to_string(ecc_len)
        );
    }

public:
    const PyGFi16 gf;
    const RSi16v<1> rs16;
    const RSi16v<4> rs16x4;
    const RSi16v<8> rs16x8;
    const RSi16v<16> rs16x16;
    bool simd_x4;
    bool simd_x8;
    bool simd_x16;

    const size_t block_len;
    const size_t message_len;
    const size_t ecc_len;
    const size_t interleave;
    const size_t interleaved_block_len;
    const size_t interleaved_message_len;
    const size_t interleaved_ecc_len;
    const size_t repair_temp_len;
    size_t vec_align;

    inline PyRSi16(
            size_t block_len,
            uint16_t ecc_len,
            size_t interleave,
            uint32_t primitive,
            std::optional<bool> simd_x4,
            std::optional<bool> simd_x8,
            std::optional<bool> simd_x16
    ):
        PyRSi16(
            rs_data(block_len, ecc_len),
            interleave,
            primitive,
            simd_x4.value_or(__builtin_cpu_supports("sse2")),
            simd_x8.value_or(__builtin_cpu_supports("avx2")),
            simd_x16.value_or(__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"))
        )
    { }

    template<typename Src, typename Dst>
    inline void encode_blocks(const Src src[], size_t full_blocks, Dst dst[]) {
        _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            auto temp = new_aligned<GFT>(block_len * SIMD_W, SIMD_W * sizeof(GFT));

            size_t block = 0;
            for (; full_blocks - block >= SIMD_W; block += SIMD_W)
                encode_block<SIMD_W>(rs, &src[block * message_len], &temp[0], &dst[block * ecc_len]);

            auto remainder = (full_blocks - block) * message_len;
            if (remainder)
                encode_block<SIMD_W>(rs, &src[block * message_len], remainder, &temp[0], &dst[block * ecc_len]);
        });
    }

    template<typename Src, typename Dst>
    inline void encode_interleaved_blocks(const Src src[], size_t full_blocks, Dst dst[]) {
        _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            // TODO: handle case when interleave is not a multiple of SIMD_W
            auto temp = new_aligned<GFT>(block_len * SIMD_W, SIMD_W * sizeof(GFT));
            for (size_t block = 0; block < full_blocks; ++block)
                _encode_interleaved<SIMD_W>(rs, &src[0], &temp[0], &dst[0]);
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

    inline void repair_block(uint32_t block[], std::vector<size_t> const& error_pos, uint32_t temp_ntt1_ecc6[]) const {
        auto error_pos_rbo = std::vector<size_t>(error_pos.size());
        for (size_t i = 0; i < error_pos.size(); ++i)
            error_pos_rbo[i] = rs16.rbo(error_pos[i]);

        rs16.repair(&block[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp_ntt1_ecc6[0]);
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

        py::class_<PyRSi16>(m, "RSi16")
            .def_property_readonly("block_len", [](PyRSi16& self) { return self.interleaved_block_len; }, R"(Interleaved block length in number of elements)")
            .def_property_readonly("block_size", [](PyRSi16& self) { return self.interleaved_block_len * sizeof(uint16_t); }, R"(Interleaved block size in bytes)")
            .def_property_readonly("message_len", [](PyRSi16& self) { return self.interleaved_message_len; }, R"(Interleaved message length in number of elements)")
            .def_property_readonly("message_size", [](PyRSi16& self) { return self.interleaved_message_len * sizeof(uint16_t); }, R"(Interleaved message size in bytes)")
            .def_property_readonly("ecc_len", [](PyRSi16& self) { return self.interleaved_ecc_len; }, R"(Interleaved error correction code length in number of elements)")
            .def_property_readonly("ecc_size", [](PyRSi16& self) { return self.interleaved_ecc_len * sizeof(uint16_t); }, R"(Interleaved error correction code size in bytes)")
            .def_property_readonly("rs_block_len", [](PyRSi16& self) { return self.block_len; }, R"(Reed-Solomon block length in number of elements)")
            .def_property_readonly("rs_block_size", [](PyRSi16& self) { return self.block_len * sizeof(uint16_t); }, R"(Reed-Solomon block size in bytes)")
            .def_property_readonly("rs_message_len", [](PyRSi16& self) { return self.message_len; }, R"(Reed-Solomon message length in number of elements)")
            .def_property_readonly("rs_message_size", [](PyRSi16& self) { return self.message_len * sizeof(uint16_t); }, R"(Reed-Solomon message size in bytes)")
            .def_property_readonly("rs_ecc_len", [](PyRSi16& self) { return self.ecc_len; }, R"(Reed-Solomon error correction code length in number of elements)")
            .def_property_readonly("rs_ecc_size", [](PyRSi16& self) { return self.ecc_len * sizeof(uint16_t); }, R"(Reed-Solomon error correction code size in bytes)")
            .def_property_readonly("ntt_len", [](PyRSi16& self) { return self.rs16.ntt_len; }, R"(Number theoretic transform length)")
            .def_property_readonly("ntt_size", [](PyRSi16& self) { return self.rs16.ntt_len * sizeof(uint16_t); }, R"(Number theoretic transform size in bytes)")
            .def_property_readonly("root", [](PyRSi16& self) { return self.rs16.root; }, R"(:math:`n`-th root of unity used by NTT)")
            .def_property_readonly("gf", [](PyRSi16& self) -> auto const& { return self.gf; }, R"(:class:`GFi16` instance for :math:`GF(65537)`)")
            .def_property_readonly("interleave", [](PyRSi16& self) { return self.interleave; }, R"(Number of interleaved columns)")

            .def_property_readonly("simd_x4", [](PyRSi16& self) { return self.simd_x4; }, R"(SIMD x4 encoding enabled (SSE2))")
            .def_property_readonly("simd_x8", [](PyRSi16& self) { return self.simd_x8; }, R"(SIMD x8 encoding enabled (AVX2))")
            .def_property_readonly("simd_x16", [](PyRSi16& self) { return self.simd_x16; }, R"(SIMD x16 encoding enabled (AVX512))")

            .def(py::init<
                    size_t,    // block_len
                    uint16_t,  // ecc_len
                    size_t,    // interleave
                    uint32_t,  // primitive
                    std::optional<bool>,
                    std::optional<bool>,
                    std::optional<bool>
                >(),
                R"(Instantiate a Reed-Solomon encoder with the given configuration)",
                "block_len"_a,
                "ecc_len"_a,
                "interleave"_a = 1,
                py::kw_only(),
                "primitive"_a = 3,
                "simd_x4"_a = py::none(),
                "simd_x8"_a = py::none(),
                "simd_x16"_a = py::none()
            )

            .def("__sizeof_cpp__", [](PyRSi16& self) { return sizeof(self); })

            .def("encode", cast_args(&PyRSi16::py_encode),
                R"(Systematic encode)",
                "buffer"_a)

            .def("synd", cast_args(&PyRSi16::py_synd),
                R"(Calculate syndromes for the given message and ecc buffers)",
                "message"_a,
                "ecc"_a)

            .def("repair", cast_args(&PyRSi16::py_repair),
                R"(Repair a block with the given error locations)",
                "message"_a,
                "ecc"_a,
                "error_pos"_a = py::none())

            .def("rbo", [](PyRSi16& self, uint16_t i) { return self.rs16.rbo(i); },
                R"(Reverse Bit Order)",
                "i"_a)

            .doc() = R"(Reed-Solomon coding over :math:`GF(65537)`)";
    }

private:
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

    template<bool Prefetch = true, size_t PrefetchDistance = 0, typename Src, typename Dst>
    inline void copy_stride(const Src *src, size_t src_stride, Dst *dst, size_t dst_stride, size_t cols, size_t count) const {
        size_t src_i = 0;
        size_t dst_i = 0;
        for (size_t i = 0; i < count; ++i) {
            if constexpr (Prefetch) {
                auto next = src_i + src_stride;
                __builtin_prefetch(&src[next]);

                if constexpr (PrefetchDistance)
                    __builtin_prefetch(&src[next + PrefetchDistance]);
            }
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
            copy_stride<true, SIMD_W / 2>(&src[i * SIMD_W], interleave, &temp[0], SIMD_W, SIMD_W, message_len);
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
            copy_stride<true, SIMD_W / 2>(&src[encoded_cols], interleave, &temp[0], SIMD_W, remaining_cols, message_len);
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

        auto temp_ntt1_ecc6 = new_aligned<GFT>((repair_temp_len) * SIMD_W, SIMD_W * sizeof(GFT));
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

            rs.repair(&buf[0], &temp_ntt1_ecc6[0]);

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

            rs.repair(&buf[0], &temp_ntt1_ecc6[0]);

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

        auto temp_ntt1_ecc6 = new_aligned<GFT>((repair_temp_len) * SIMD_W, SIMD_W * sizeof(GFT));
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

            rs.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp_ntt1_ecc6[0]);

            copy_stride(&buf[0], SIMD_W, &message[i * SIMD_W], interleave, SIMD_W, message_len);

            // Sequential ecc
            // copy_transposed(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W * ecc_len], ecc_len);

            // Interleaved ecc
            copy_stride(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W], interleave, SIMD_W, ecc_len);
        }

        size_t encoded_cols = vec_cols * SIMD_W;
        if (encoded_cols < col_count) {
            size_t remaining_cols = col_count - encoded_cols;
            // TODO: fill_stride
            std::fill_n(&buf[0], block_len * SIMD_W, GFT{0});
            copy_stride(&message[encoded_cols], interleave, &buf[0], SIMD_W, remaining_cols, message_len);

            // Sequential ecc
            // copy_transposed(&ecc[encoded_cols * ecc_len], ecc_len, ecc_len, &buf[message_len * SIMD_W], SIMD_W, remaining_cols);

            // Interleaved ecc
            copy_stride(&ecc[encoded_cols], interleave, &buf[message_len * SIMD_W], SIMD_W, remaining_cols, ecc_len);

            rs.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp_ntt1_ecc6[0]);

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

        py_assert(buf.size % interleaved_message_len == 0, std::to_string(buf.size));

        size_t full_blocks = buf.size / interleaved_message_len;
        size_t output_size = full_blocks * interleaved_ecc_len;


        auto output = py::bytearray(nullptr, output_size * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

        if (interleave == 1)
            encode_blocks(&buf[0], full_blocks, &output_data[0]);
        else
            encode_interleaved_blocks(&buf[0], full_blocks, &output_data[0]);

        return output;
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

            auto temp_ntt1_ecc6 = new_aligned<GFT>(repair_temp_len, sizeof(GFT));
            if (error_pos) {
                py_assert(error_pos->size() <= ecc_len);

                auto error_pos_rbo = std::vector<size_t>(error_pos->size());
                for (size_t i = 0; i < error_pos->size(); ++i)
                    error_pos_rbo[i] = rs16.rbo((*error_pos)[i]);

                rs16.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp_ntt1_ecc6[0]);

            } else {
                rs16.repair(&buf[0], &temp_ntt1_ecc6[0]);
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
