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
#include "rsi16v.hpp"
#include "pyntt.hpp"

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
        GFT primitive,
        bool simd_x4,
        bool simd_x8,
        bool simd_x16
    ):
        gf(primitive),
        ntt(gf, args.block_len, args.ecc_len),
        rs16(ntt, args.block_len, args.ecc_len),
        rs16x4(ntt, args.block_len, args.ecc_len),
        rs16x8(ntt, args.block_len, args.ecc_len),
        rs16x16(ntt, args.block_len, args.ecc_len),
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
        repair_temp_len(ntt.ntt_len + ecc_len * 6)
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
    const NTT ntt;
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
            GFT primitive,
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
    inline void encode_blocks(const Src src[], size_t full_blocks, Dst dst[]) const {
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
    inline void encode_interleaved_blocks(const Src src[], size_t full_blocks, Dst dst[]) const {
        _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            // TODO: handle case when interleave is not a multiple of SIMD_W
            auto temp = new_aligned<GFT>(block_len * SIMD_W, SIMD_W * sizeof(GFT));
            for (size_t block = 0; block < full_blocks; ++block)
                _encode_interleaved<SIMD_W>(rs, &src[0], &temp[0], &dst[0]);
        });
    }

    template<typename Src, typename Dst>
    inline void encode_interleaved(const Src src[], Dst dst[]) const {
        _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            auto temp = new_aligned<GFT>(block_len * SIMD_W, SIMD_W * sizeof(GFT));
            _encode_interleaved<SIMD_W>(rs, &src[0], &temp[0], &dst[0]);
        });
    }

    template<typename Msg, typename Ecc>
    inline RepairStatus repair_interleaved(Msg message[], Ecc ecc[], size_t col_start, size_t col_count) const {
        return _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            return _repair_interleaved<SIMD_W>(rs, &message[0], &ecc[0], col_start, col_count);
        });
    }

    template<typename Msg, typename Ecc>
    inline RepairStatus repair_interleaved(Msg message[], Ecc ecc[], size_t col_start, size_t col_count, std::vector<size_t> const& error_pos) const {
        auto error_pos_rbo = std::vector<size_t>(error_pos.size());
        for (size_t i = 0; i < error_pos.size(); ++i)
            error_pos_rbo[i] = ntt.rbo(error_pos[i]);

        return _simd_dispatch([&]<size_t SIMD_W>(std::integral_constant<size_t, SIMD_W>, auto& rs) {
            return _repair_interleaved<SIMD_W>(rs, &message[0], &ecc[0], col_start, col_count, error_pos_rbo);
        });
    }

    inline void repair_block(GFT block[], std::vector<size_t> const& error_pos, GFT temp_ntt1_ecc6[]) const {
        auto error_pos_rbo = std::vector<size_t>(error_pos.size());
        for (size_t i = 0; i < error_pos.size(); ++i)
            error_pos_rbo[i] = ntt.rbo(error_pos[i]);

        rs16.repair(&block[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp_ntt1_ecc6[0]);
    }

    inline void synd_block(GFT block[]) const {
        ntt.pntt(&block[0]);
    }

    template<typename Msg, typename Ecc>
    inline void synd_blocks(const Msg msg[], const Ecc ecc[], size_t count, GFT temp[], GFT synds[]) const {
        // len(temp) == block_len
        // len(synds) == count * ecc_len
        for (size_t i = 0; i < count; ++i) {
            std::copy_n(&msg[i * message_len], message_len, &temp[0]);
            std::copy_n(&ecc[i * ecc_len], ecc_len, &temp[message_len]);
            ntt.pntt(&temp[0]);
            std::copy_n(&temp[0], ecc_len, &synds[i * ecc_len]);
        }
    }

    inline void synd_interleaved(const uint16_t msg[], const uint16_t ecc[], GFT temp[], GFT synds[]) {
        // len(temp) == block_len
        // len(synds) == count * ecc_len

        const size_t SIMD_W = 1;
        for (size_t i = 0; i < interleave; ++i) {
            // std::copy_n(&msg[i * message_len], message_len, &temp[0]);
            vec::copy_stride(&msg[i * SIMD_W], interleave, &temp[0], SIMD_W, SIMD_W, message_len);
            // std::copy_n(&ecc[i * ecc_len], ecc_len, &temp[message_len]);
            vec::copy_transposed(&ecc[i * ecc_len], ecc_len, &temp[message_len * SIMD_W], SIMD_W);
            ntt.pntt(&temp[0]);
            // std::copy_n(&temp[0], ecc_len, &synds[i * ecc_len]);
            vec::copy_transposed(&temp[0], SIMD_W, &synds[i * ecc_len * SIMD_W], ecc_len);
        }
    }

    inline void mix_ecc(GFT ecc[]) const {
        rs16.mix_ecc(&ecc[0]);
    }

    inline size_t message_offset(size_t row, size_t col = 0) {
        py_assert(row < message_len);
        py_assert(col < interleave);
        return row * interleave + col;
    }

    inline size_t ecc_offset(size_t row, size_t col = 0) {
        py_assert(row < ecc_len);
        py_assert(col < interleave);
        return row * interleave + col;
    }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::enum_<RepairStatus>(m, "RepairStatus", R"(Repair operation status)")
            .value("NoErrors", RepairStatus::NoErrors)
            .value("NoErrorsZero", RepairStatus::NoErrorsZero)
            .value("RepairOk", RepairStatus::RepairOk)
            .value("RepairFail", RepairStatus::RepairFail)
            .value("ErrorLocationFail", RepairStatus::ErrorLocationFail)
        ;

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
            .def_property_readonly("ntt_len", [](PyRSi16& self) { return self.ntt.ntt_len; }, R"(Number theoretic transform length)")
            .def_property_readonly("ntt_size", [](PyRSi16& self) { return self.ntt.ntt_len * sizeof(uint16_t); }, R"(Number theoretic transform size in bytes)")
            .def_property_readonly("root", [](PyRSi16& self) { return self.rs16.root; }, R"(:math:`n`-th root of unity used by NTT)")
            .def_property_readonly("gf", [](PyRSi16& self) -> auto const * { return &self.gf; }, R"(:class:`GFi16` instance for :math:`GF(65537)`)")
            .def_property_readonly("ntt", [](PyRSi16& self) -> auto const * { return &self.ntt; }, R"(:class:`NTT` instance)")
            .def_property_readonly("interleave", [](PyRSi16& self) { return self.interleave; }, R"(Number of interleaved columns)")

            .def_property_readonly("simd_x4", [](PyRSi16& self) { return self.simd_x4; }, R"(SIMD x4 encoding enabled (SSE2))")
            .def_property_readonly("simd_x8", [](PyRSi16& self) { return self.simd_x8; }, R"(SIMD x8 encoding enabled (AVX2))")
            .def_property_readonly("simd_x16", [](PyRSi16& self) { return self.simd_x16; }, R"(SIMD x16 encoding enabled (AVX512))")

            .def(py::init<
                    size_t,    // block_len
                    uint16_t,  // ecc_len
                    size_t,    // interleave
                    GFT,  // primitive
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

            .def("repair", cast_args(&PyRSi16::py_repair),
                R"(Repair a block with the given error locations)",
                "message"_a,
                "ecc"_a,
                "error_pos"_a = py::none())

            .def("_synd", cast_args(&PyRSi16::py_synd),
                R"(Calculate syndromes for the given message and ecc buffers)",
                "message"_a,
                "ecc"_a)

            .def("_sugiyama", &PyRSi16::py_sugiyama<1>, R"(Compute error locator and evaluator polynomials)", "synd"_a)
            .def("_sugiyama4", &PyRSi16::py_sugiyama<4>, R"(Compute error locator and evaluator polynomials)", "synd"_a)
            .def("_sugiyama8", &PyRSi16::py_sugiyama<8>, R"(Compute error locator and evaluator polynomials)", "synd"_a)
            .def("_sugiyama16", &PyRSi16::py_sugiyama<16>, R"(Compute error locator and evaluator polynomials)", "synd"_a)

            .def("message_offset", &PyRSi16::message_offset, R"(Calculate message offset in number of elements)", "row"_a, "col"_a)
            .def("ecc_offset", &PyRSi16::ecc_offset, R"(Calculate ECC offset in number of elements)", "row"_a, "col"_a)

            .def("_roots", &PyRSi16::py_roots,
                R"(Find locator polynomial roots)",
                "synd"_a)

            .doc() = R"(Reed-Solomon coding over :math:`GF(65537)`)"
        ;
    }

private:
    template<size_t SIMD_W, typename RS, typename Src, typename Dst>
    inline void _encode_interleaved(RS const& rs, const Src src[], GFT temp[], Dst dst[]) const {
        // src_size = message_len * interleave
        // dst_size = ecc_len * interleave
        // block_len = message_len + ecc_len
        // chunk_size = block_len * interleave
        // temp_size = block_len * SIMD_W
        size_t vec_cols = interleave / SIMD_W;
        for (size_t i = 0; i < vec_cols; ++i) {
            vec::copy_stride<true, SIMD_W / 2>(&src[i * SIMD_W], interleave, &temp[0], SIMD_W, SIMD_W, message_len);
            std::fill_n(&temp[message_len * SIMD_W], ecc_len * SIMD_W, 0);
            // std::memset(&temp[message_len * SIMD_W], 0, ecc_len * SIMD_W * sizeof(GFT));

            rs.encode(&temp[0]);
            // std::copy_n(&temp[0], ecc_len * SIMD_W, &dst[i * SIMD_W * ecc_len]);

            // Sequential ecc
            // vec::copy_transposed(&temp[0], SIMD_W, &dst[i * ecc_len * SIMD_W], ecc_len);

            // Interleaved ecc
            vec::copy_stride(&temp[0], SIMD_W, &dst[i * SIMD_W], interleave, SIMD_W, ecc_len);
        }
        size_t encoded_cols = vec_cols * SIMD_W;
        if (encoded_cols < interleave) {
            size_t remaining_cols = interleave - encoded_cols;
            vec::copy_stride<true, SIMD_W / 2>(&src[encoded_cols], interleave, &temp[0], SIMD_W, remaining_cols, message_len);
            std::fill_n(&temp[message_len * SIMD_W], ecc_len * SIMD_W, 0);
            // std::memset(&temp[message_len * SIMD_W], 0, ecc_len * SIMD_W * sizeof(GFT));

            rs.encode(&temp[0]);
            // vec::copy_stride(&temp[0], SIMD_W, &dst[encoded_cols], remaining_cols, remaining_cols, ecc_len);

            // Sequential ecc
            // vec::copy_transposed(&temp[0], SIMD_W, remaining_cols, &dst[encoded_cols * ecc_len], ecc_len, ecc_len);

            // Interleaved ecc
            vec::copy_stride(&temp[0], SIMD_W, &dst[encoded_cols], interleave, remaining_cols, ecc_len);
        }
    }

    template<size_t SIMD_W, typename T, typename Src, typename Dst>
    inline void encode_block(T const& rs, const Src src[], GFT temp[], Dst dst[]) const {
        vec::copy_transposed(&src[0], message_len, &temp[0], SIMD_W);
        std::fill_n(&temp[message_len * SIMD_W], ecc_len * SIMD_W, 0);
        // std::memset(&temp[message_len * SIMD_W], 0, ecc_len * SIMD_W * sizeof(GFT));

        rs.encode(&temp[0]);
        vec::copy_transposed(&temp[0], SIMD_W, &dst[0], ecc_len);
    }

    template<size_t SIMD_W, typename T, typename Src, typename Dst>
    inline void encode_block(T const& rs, const Src src[], size_t src_size, GFT temp[], Dst dst[]) const {
        std::fill_n(&temp[0], block_len * SIMD_W, 0);

        auto cols = src_size / message_len;
        vec::copy_transposed(&src[0], message_len, message_len, &temp[0], SIMD_W, cols);

        size_t remainder = src_size - cols * message_len;
        if (remainder) {
            vec::copy_transposed(&src[src_size - remainder], remainder, remainder, &temp[cols], SIMD_W, 1);
            ++cols;
        }

        rs.encode(&temp[0]);
        vec::copy_transposed(&temp[0], SIMD_W, cols, &dst[0], ecc_len, ecc_len);
    }

    template<size_t SIMD_W, typename Rs, typename Msg, typename Ecc>
    inline RepairStatus _repair_interleaved(Rs const& rs, Msg message[], Ecc ecc[], size_t col_start, size_t col_count) const {
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

        RepairStatus res = RepairStatus::NoErrors;

        size_t vec_cols = col_count / SIMD_W;
        for (size_t i = 0; i < vec_cols; ++i) {
            vec::copy_stride(&message[i * SIMD_W], interleave, &buf[0], SIMD_W, SIMD_W, message_len);

            // Sequential ecc
            // vec::copy_transposed(&ecc[i * SIMD_W * ecc_len], ecc_len, &buf[message_len * SIMD_W], SIMD_W);

            // Interleaved ecc
            vec::copy_stride(&ecc[i * SIMD_W], interleave, &buf[message_len * SIMD_W], SIMD_W, SIMD_W, ecc_len);

            res = std::max(res, rs.repair(&buf[0], &temp_ntt1_ecc6[0]));

            vec::copy_stride(&buf[0], SIMD_W, &message[i * SIMD_W], interleave, SIMD_W, message_len);

            // Sequential ecc
            // vec::copy_transposed(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W * ecc_len], ecc_len);

            // Interleaved ecc
            vec::copy_stride(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W], interleave, SIMD_W, ecc_len);
        }

        size_t encoded_cols = vec_cols * SIMD_W;
        if (encoded_cols < col_count) {
            std::fill_n(&buf[0], block_len * SIMD_W, GFT{0});

            size_t remaining_cols = col_count - encoded_cols;
            vec::copy_stride(&message[encoded_cols], interleave, &buf[0], SIMD_W, remaining_cols, message_len);

            // Sequential ecc
            // vec::copy_transposed(&ecc[encoded_cols * ecc_len], ecc_len, ecc_len, &buf[message_len * SIMD_W], SIMD_W, remaining_cols);

            // Interleaved ecc
            vec::copy_stride(&ecc[encoded_cols], interleave, &buf[message_len * SIMD_W], SIMD_W, remaining_cols, ecc_len);

            res = std::max(res, rs.repair(&buf[0], &temp_ntt1_ecc6[0]));

            vec::copy_stride(&buf[0], SIMD_W, &message[encoded_cols], interleave, remaining_cols, message_len);

            // Sequential ecc
            // vec::copy_transposed(&buf[message_len * SIMD_W], SIMD_W, remaining_cols, &ecc[encoded_cols * ecc_len], ecc_len, ecc_len);

            // Interleaved ecc
            vec::copy_stride(&buf[message_len * SIMD_W], SIMD_W, &ecc[encoded_cols], interleave, remaining_cols, ecc_len);
        }
        return res;
    }

    template<size_t SIMD_W, typename Rs, typename Msg, typename Ecc>
    inline RepairStatus _repair_interleaved(Rs const& rs, Msg message[], Ecc ecc[], size_t col_start, size_t col_count, std::vector<size_t> const& error_pos_rbo) const {
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

        RepairStatus res = RepairStatus::NoErrors;

        size_t vec_cols = col_count / SIMD_W;
        for (size_t i = 0; i < vec_cols; ++i) {
            vec::copy_stride(&message[i * SIMD_W], interleave, &buf[0], SIMD_W, SIMD_W, message_len);

            // Sequential ecc
            // vec::copy_transposed(&ecc[i * SIMD_W * ecc_len], ecc_len, &buf[message_len * SIMD_W], SIMD_W);

            // Interleaved ecc
            vec::copy_stride(&ecc[i * SIMD_W], interleave, &buf[message_len * SIMD_W], SIMD_W, SIMD_W, ecc_len);

            res = std::max(res, rs.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp_ntt1_ecc6[0]));

            vec::copy_stride(&buf[0], SIMD_W, &message[i * SIMD_W], interleave, SIMD_W, message_len);

            // Sequential ecc
            // vec::copy_transposed(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W * ecc_len], ecc_len);

            // Interleaved ecc
            vec::copy_stride(&buf[message_len * SIMD_W], SIMD_W, &ecc[i * SIMD_W], interleave, SIMD_W, ecc_len);
        }

        size_t encoded_cols = vec_cols * SIMD_W;
        if (encoded_cols < col_count) {
            size_t remaining_cols = col_count - encoded_cols;
            // TODO: fill_stride
            std::fill_n(&buf[0], block_len * SIMD_W, GFT{0});
            vec::copy_stride(&message[encoded_cols], interleave, &buf[0], SIMD_W, remaining_cols, message_len);

            // Sequential ecc
            // vec::copy_transposed(&ecc[encoded_cols * ecc_len], ecc_len, ecc_len, &buf[message_len * SIMD_W], SIMD_W, remaining_cols);

            // Interleaved ecc
            vec::copy_stride(&ecc[encoded_cols], interleave, &buf[message_len * SIMD_W], SIMD_W, remaining_cols, ecc_len);

            res = std::max(res, rs.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp_ntt1_ecc6[0]));

            vec::copy_stride(&buf[0], SIMD_W, &message[encoded_cols], interleave, remaining_cols, message_len);

            // Sequential ecc
            // vec::copy_transposed(&buf[message_len * SIMD_W], SIMD_W, remaining_cols, &ecc[encoded_cols * ecc_len], ecc_len, ecc_len);

            // Interleaved ecc
            vec::copy_stride(&buf[message_len * SIMD_W], SIMD_W, &ecc[encoded_cols], interleave, remaining_cols, ecc_len);
        }
        return res;
    }

    template<typename F>
    inline auto _simd_dispatch(F&& f) const {
        if (simd_x16)
            return f(std::integral_constant<size_t, 16>{}, rs16x16);
        else if (simd_x8)
            return f(std::integral_constant<size_t, 8>{}, rs16x8);
        else if (simd_x4)
            return f(std::integral_constant<size_t, 4>{}, rs16x4);
        else
            return f(std::integral_constant<size_t, 1>{}, rs16);
    }

    template<size_t W>
    inline RSi16v<W> const& rs() const {
        if constexpr (W == 16)
            return rs16x16;
        else if constexpr (W == 8)
            return rs16x8;
        else if constexpr (W == 4)
            return rs16x4;
        else if constexpr (W == 1)
            return rs16;
        else
            static_assert(W == 1 || W == 4 || W == 8 || W == 16, "Unsupported SIMD width");
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

    inline RepairStatus py_repair(buffer_rw<uint16_t> message, buffer_rw<uint16_t> ecc, std::optional<std::vector<size_t>> const& error_pos) {
        RepairStatus res = RepairStatus::NoErrors;

        if (error_pos && error_pos->empty())
            return RepairStatus::NoErrors;

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
                    error_pos_rbo[i] = ntt.rbo((*error_pos)[i]);

                res = rs16.repair(&buf[0], &error_pos_rbo[0], error_pos_rbo.size(), &temp_ntt1_ecc6[0]);
            } else {
                res = rs16.repair(&buf[0], &temp_ntt1_ecc6[0]);
            }

            std::copy_n(&buf[0], message_len, &message[0]);
            std::copy_n(&buf[message_len], ecc_len, &ecc[0]);
        } else {
            py_assert(message.size == interleaved_message_len, std::to_string(message.size));
            py_assert(ecc.size == interleaved_ecc_len, std::to_string(ecc.size));

            if (error_pos) {
                py_assert(error_pos->size() <= ecc_len);
                res = repair_interleaved(&message[0], &ecc[0], 0, interleave, *error_pos);
            } else {
                res = repair_interleaved(&message[0], &ecc[0], 0, interleave);
            }
        }

        return res;
    }

    inline std::vector<GFT> py_synd(buffer_ro<uint16_t> message, buffer_ro<uint16_t> ecc) {
        py_assert(message.size % message_len == 0, std::to_string(message.size));
        py_assert(ecc.size % ecc_len == 0, std::to_string(ecc.size));

        size_t count = message.size / message_len;
        py_assert(ecc.size / ecc_len == count);

        auto temp = new_aligned<GFT>(block_len, sizeof(GFT));
        std::vector<GFT> synd;
        synd.resize(count * ecc_len);

        std::copy_n(&message[0], message_len, &temp[0]);
        std::copy_n(&ecc[0], ecc_len, &temp[message_len]);

        ntt.pntt(&temp[0]);
        synd_blocks(&message[0], &ecc[0], count, &temp[0], &synd[0]);

        return synd;
    }

    template<size_t W>
    inline py::tuple py_sugiyama(std::vector<GFT> synd) {
        auto repair_temp = new_aligned<GFT>(repair_temp_len * W, W * sizeof(GFT));
        auto const& rs = this->rs<W>();

        auto a1 = &repair_temp[0];
        auto r1 = &repair_temp[W * ecc_len];
        auto temp_ecc4 = &repair_temp[2 * W * ecc_len];

        vec::copy_transposed(&synd[0], ecc_len, &r1[0], W);

        rs.sugiyama(&a1[0], &r1[0], &temp_ecc4[0]);

        std::vector<GFT> locator, evaluator;
        locator.resize(W * ecc_len);
        evaluator.resize(W * ecc_len);

        vec::copy_transposed(&a1[0], W, &locator[0], ecc_len);
        vec::copy_transposed(&r1[0], W, &evaluator[0], ecc_len);

        return py::make_tuple(locator, evaluator);
    }

    inline std::vector<size_t> py_roots(std::vector<GFT> poly) {
        py_assert(poly.size() == ecc_len);
        auto temp = new_aligned<GFT>(block_len, sizeof(GFT));

        std::copy_n(&poly[0], ecc_len, &temp[0]);
        std::fill_n(&temp[ecc_len], block_len - ecc_len, 0);

        ntt.pintt_ecc(&temp[0]);

        std::vector<size_t> zeros;
        for (size_t i = 0; i < block_len; ++i) {
            if (temp[i] == 0)
                zeros.push_back(i);
        }
        return zeros;
    }
};
