/**************************************************************************
 * pyCIRC16.hpp
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

#pragma once

#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pylogging.hpp"
#include "pyrsi16.hpp"


namespace py = pybind11;


class PyCIRC16 {
private:
    PyRSi16 rsi;
    PyRSi16 rso;

public:
    const size_t block_len;
    const size_t message_len;
    const size_t rsi_interleaved_ecc_len;
    const size_t rsio_ecc_len;
    const size_t ecc_len;
    const size_t interleave;

    inline PyCIRC16(
        size_t inner_block_len,
        size_t inner_ecc_len,
        size_t outer_block_len,
        size_t outer_ecc_len,
        size_t interleave,
        size_t primitive,
        std::optional<bool> simd_x4,
        std::optional<bool> simd_x8,
        std::optional<bool> simd_x16
    ):
        rsi(inner_block_len, inner_ecc_len, 1, primitive, simd_x4, simd_x8, simd_x16),
        rso(outer_block_len, outer_ecc_len, rsi.message_len * interleave, primitive, simd_x4, simd_x8, simd_x16),
        block_len(inner_block_len * outer_block_len * interleave),
        message_len(rso.interleaved_message_len),
        rsi_interleaved_ecc_len(rsi.ecc_len * rso.message_len * interleave),
        rsio_ecc_len(rsi.ecc_len * rso.ecc_len * interleave),
        ecc_len(rso.interleaved_ecc_len + rsi_interleaved_ecc_len + rsio_ecc_len),
        interleave(interleave)
    { }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyCIRC16>(m, "CIRC16")
            .def(py::init<
                    size_t,  // inner_block_len
                    size_t,  // inner_ecc_len
                    size_t,  // outer_block_len
                    size_t,  // outer_ecc_len
                    size_t,  // interleave
                    GFT,  // primitive
                    std::optional<bool>,  // simd_x4
                    std::optional<bool>,  // simd_x8
                    std::optional<bool>   // simd_x16
                >(),
                R"(Cross-interleaved Reed-Solomon coder)",
                "inner_block_len"_a,
                "inner_ecc_len"_a,
                "outer_block_len"_a,
                "outer_ecc_len"_a,
                "interleave"_a = 1,
                py::kw_only(),
                "primitive"_a = 3,
                "simd_x4"_a = py::none(),
                "simd_x8"_a = py::none(),
                "simd_x16"_a = py::none()
            )

            .def_property_readonly("rsi", [](PyCIRC16& self) -> auto const& { return self.rsi; }, R"(Inner Reed-Solomon codec (:class:`RSi16`))")
            .def_property_readonly("rso", [](PyCIRC16& self) -> auto const& { return self.rso; }, R"(Outer Reed-Solomon codec (:class:`RSi16`))")

            .def_property_readonly("block_len", [](PyCIRC16& self) { return self.block_len; }, R"(Total block length in number of elements)")
            .def_property_readonly("block_size", [](PyCIRC16& self) { return self.block_len * sizeof(uint16_t); }, R"(Total block size in bytes)")
            .def_property_readonly("message_len", [](PyCIRC16& self) { return self.message_len; }, R"(Total message length in number of elements)")
            .def_property_readonly("message_size", [](PyCIRC16& self) { return self.message_len * sizeof(uint16_t); }, R"(Total message size in bytes)")
            .def_property_readonly("ecc_len", [](PyCIRC16& self) { return self.ecc_len; }, R"(Total ECC length in number of elements)")
            .def_property_readonly("ecc_size", [](PyCIRC16& self) { return self.ecc_len * sizeof(uint16_t); }, R"(Total ECC size in bytes)")
            .def_property_readonly("inner_block_len", [](PyCIRC16& self) { return self.rsi.block_len; }, R"(Inner codec block length in number of elements)")
            .def_property_readonly("inner_block_size", [](PyCIRC16& self) { return self.rsi.block_len * sizeof(uint16_t); }, R"(Inner codec block size in bytes)")
            .def_property_readonly("inner_message_len", [](PyCIRC16& self) { return self.rsi.message_len; }, R"(Inner codec message length in number of elements)")
            .def_property_readonly("inner_message_size", [](PyCIRC16& self) { return self.rsi.message_len * sizeof(uint16_t); }, R"(Inner codec message size in bytes)")
            .def_property_readonly("inner_ecc_len", [](PyCIRC16& self) { return self.rsi.ecc_len; }, R"(Inner codec ECC length in number of elements)")
            .def_property_readonly("inner_ecc_size", [](PyCIRC16& self) { return self.rsi.ecc_len * sizeof(uint16_t); }, R"(Inner codec ECC size in bytes)")
            .def_property_readonly("outer_block_len", [](PyCIRC16& self) { return self.rso.block_len; }, R"(Outer codec block length in number of elements)")
            .def_property_readonly("outer_message_len", [](PyCIRC16& self) { return self.rso.message_len; }, R"(Outer codec message length in number of elements)")
            .def_property_readonly("outer_ecc_len", [](PyCIRC16& self) { return self.rso.ecc_len; }, R"(Outer codec ECC length in number of elements)")
            .def_property_readonly("outer_interleave", [](PyCIRC16& self) { return self.rso.interleave; }, R"(Outer codec interleaving (``inner_message_len * interleave``))")
            .def_property_readonly("interleave", [](PyCIRC16& self) { return self.interleave; }, R"(Number of interleaved CIRC blocks)")
            .def_property_readonly("primitive", [](PyCIRC16& self) { return self.rsi.gf.primitive; }, R"(Primitive value used to generate the Galois field)")
            .def("__sizeof_cpp__", [](PyCIRC16& self) { return sizeof(self); }, R"(Size of object in bytes)")
            .def_property_readonly("simd_x4", [](PyCIRC16& self) { return self.rsi.simd_x4; }, R"(SIMD x4 encoding enabled (SSE2))")
            .def_property_readonly("simd_x8", [](PyCIRC16& self) { return self.rsi.simd_x8; }, R"(SIMD x8 encoding enabled (AVX2))")
            .def_property_readonly("simd_x16", [](PyCIRC16& self) { return self.rsi.simd_x16; }, R"(SIMD x16 encoding enabled (AVX512))")

            .def("encode", cast_args(&PyCIRC16::py_encode), R"(Encode data)", "buffer"_a)
            .def("repair", cast_args(&PyCIRC16::py_repair), R"(Repair data)", "message"_a, "ecc"_a)
            .doc() = R"(Cross-interleaved Reed-Solomon coding over :math:`GF(65537)`)";
            ;
    }

private:
    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        py_assert(buf.size % message_len == 0, std::to_string(buf.size));

        size_t full_blocks = buf.size / message_len;

        auto output = py::bytearray(nullptr, full_blocks * ecc_len * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));
        auto temp = new_aligned<GFT>(rso.interleaved_ecc_len, rsi.vec_align);

        for (size_t i = 0; i < full_blocks; ++i)
            encode_block(&buf[i * message_len], &output_data[i * ecc_len]);

        return output;
    }

    inline void encode_block(const uint16_t src[], uint16_t dst[]) {
        auto temp = new_aligned<GFT>(rso.interleaved_ecc_len, rsi.vec_align);

        auto rso_ecc = &dst[0];  // size = rso.interleaved_ecc_len
        auto rsi_ecc = &dst[rso.interleaved_ecc_len];  // size = rsi_interleaved_ecc_len
        auto rsio_ecc = &dst[rso.interleaved_ecc_len + rsi_interleaved_ecc_len];  // size = rsio_ecc_len

        rso.encode_interleaved(&src[0], &temp[0]);
        std::copy_n(&temp[0], rso.interleaved_ecc_len, &rso_ecc[0]);
        rsi.encode_blocks(&src[0], rso.message_len * interleave, &rsi_ecc[0]);
        rsi.encode_blocks(&temp[0], rso.ecc_len * interleave, &rsio_ecc[0]);
    }

    inline bool py_repair(buffer_rw<uint16_t> message, buffer_rw<uint16_t> ecc) {
        log_debug("message size: %s", message.size);
        log_debug("ecc size: %s", ecc.size);
        py_assert(message.size == message_len, std::to_string(message.size) + " != " + std::to_string(message_len));
        py_assert(ecc.size == ecc_len, std::to_string(ecc.size));

        size_t inner_blocks = message.size / rsi.message_len;
        py_assert(inner_blocks == rso.message_len * interleave);
        log_info("rsi blocks: %s", inner_blocks);

        auto buf_rsi = new_aligned<GFT>(rsi.block_len, rsi.vec_align);
        auto repair_temp_rsi = new_aligned<GFT>(rsi.repair_temp_len, rsi.vec_align);
        auto synds_rsi = new_aligned<GFT>(rso.block_len * rsi.ecc_len * interleave, rsi.vec_align);

        // Widen rso ecc to u32
        auto rso_ecc = new_aligned<GFT>(rso.interleaved_ecc_len, rsi.vec_align);
        std::copy_n(&ecc[0], rso.interleaved_ecc_len, &rso_ecc[0]);

        auto rsi_ecc = &ecc[rso.interleaved_ecc_len];
        auto rsio_ecc = &ecc[rso.interleaved_ecc_len + rsi_interleaved_ecc_len];

        rsi.synd_blocks(
            &message[0],
            &rsi_ecc[0],
            rso.message_len * interleave,
            &repair_temp_rsi[0],
            &synds_rsi[0]
        );
        rsi.synd_blocks(
            &rso_ecc[0],
            &rsio_ecc[0],
            rso.ecc_len * interleave,
            &repair_temp_rsi[0],
            &synds_rsi[rsi_interleaved_ecc_len]
        );

        // Repair zeroes in outer ecc
        for (size_t k = 0; k < interleave; ++k) {
            for (size_t i = 0; i < rso.ecc_len; ++i) {
                std::vector<size_t> inner_zero_locations;
                size_t rso_offset = (k + i * interleave) * rsi.message_len;
                size_t synd_offset = (k + (rso.message_len + i) * interleave) * rsi.ecc_len;

                if (std::all_of(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len], [](auto v) { return v == 0; }))
                    continue;

                for (size_t j = 0; j < rsi.message_len; ++j) {
                    if (rso_ecc[rso_offset + j] == 0) {
                        inner_zero_locations.push_back(j);
                    }
                }

                size_t rsio_offset = (k + i * interleave) * rsi.ecc_len;
                for (size_t j = 0; j < rsi.ecc_len; ++j) {
                    if (rsio_ecc[rsio_offset + j] == 0) {
                        inner_zero_locations.push_back(rsi.message_len + j);
                    }
                }

                if (inner_zero_locations.empty())
                    continue;

                if (inner_zero_locations.size() > rsi.ecc_len) {
                    log_warning("too many zeros in outer ecc: interleave:%d row:%d count:%d", k, rso.message_len + i, inner_zero_locations.size());
                    continue;
                }

                log_debug("outer ecc zero: interleave:%d row:%d locations: %s", k, rso.message_len + i, inner_zero_locations);
                log_debug(" synd: %s", std::vector(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len]));

                std::copy_n(&rso_ecc[rso_offset], rsi.message_len, &buf_rsi[0]);
                std::copy_n(&rsio_ecc[rsio_offset], rsi.ecc_len, &buf_rsi[rsi.message_len]);
                rsi.repair_block(&buf_rsi[0], inner_zero_locations, &repair_temp_rsi[0]);
                std::copy_n(&buf_rsi[0], rsi.message_len, &rso_ecc[rso_offset]);

                if (std::all_of(inner_zero_locations.begin(), inner_zero_locations.end(),
                    [&](size_t zero_pos) { return (buf_rsi[zero_pos] & 0xffff) == 0; }))
                {
                    std::fill_n(&synds_rsi[synd_offset], rsi.ecc_len, 0);
                } else {
                    log_warning("could not repair zeros in outer ecc");
                    log_warning(" interleave:%d row:%d", k, i);
                    log_warning(" locations: %s", inner_zero_locations);
                    log_warning(" synd: %s", std::vector(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len]));
                }
            }
        }

        std::vector<size_t> outer_error_locations;
        for (size_t k = 0; k < interleave; ++k) {
            for (size_t i = 0; i < rso.block_len; ++i) {
                // size_t synd_offset = k * rsi.ecc_len + i * rsi.ecc_len * interleave;
                size_t synd_offset = (k + i * interleave) * rsi.ecc_len;
                if (std::any_of(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len], [](auto v) { return v != 0; })) {
                    bool inner_ecc_has_zeros = false;
                    // Repair zeroes in rsi ecc
                    if (std::any_of(&rsi_ecc[synd_offset], &rsi_ecc[synd_offset + rsi.ecc_len], [](auto v) { return v == 0; })) {
                        inner_ecc_has_zeros = true;
                        std::copy_n(&synds_rsi[synd_offset], rsi.ecc_len, &repair_temp_rsi[0]);
                        rsi.mix_ecc(&repair_temp_rsi[0]);

                        if (std::all_of(&repair_temp_rsi[0], &repair_temp_rsi[rsi.ecc_len], [](auto v) { return (v & 0xffff) == 0; })) {
                            log_info("inner ecc zero: interleave:%d row:%d", k, i);
                            continue;
                        } else {
                            log_warning("could not repair zeros in inner ecc");
                            log_warning(" interleave:%d row:%d", k, i);
                            log_warning(" synd: %s", std::vector(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len]));
                        }
                    }

                    log_info("inner check fail: interleave:%d row:%d", k, i);
                    log_info(" inner ecc has zeros: %s", inner_ecc_has_zeros);
                    if (i > rso.message_len) {
                        size_t rso_ecc_offset = (k + (i - rso.message_len) * interleave) * rsi.message_len;
                        bool outer_ecc_has_zeros = std::any_of(&rso_ecc[rso_ecc_offset], &rso_ecc[rso_ecc_offset + rsi.message_len], [](auto v) { return v == 0; });
                        log_info(" outer ecc has zeros: %s", outer_ecc_has_zeros);

                        if (outer_ecc_has_zeros) {
                            log_info(" outer ecc: %s", std::vector(&rso_ecc[rso_ecc_offset], &rso_ecc[rso_ecc_offset + rsi.message_len]));
                            log_info(" synd: %s", std::vector(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len]));
                        }
                    }
                    log_debug(" ecc: %s", std::vector(&rsi_ecc[synd_offset], &rsi_ecc[synd_offset + rsi.ecc_len]));
                    log_debug(" synd: %s", std::vector(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len]));
                    outer_error_locations.push_back(i);
                }
            }

            if (outer_error_locations.empty())
                continue;

            log_warning("errors found: interleave:%d count:%d", k, outer_error_locations.size());
            log_debug("error locations: %s", outer_error_locations);

            if (outer_error_locations.size() <= rso.ecc_len) {
                rso.repair_interleaved(&message[0], &rso_ecc[0], k * rsi.message_len, rsi.message_len, outer_error_locations);
            } else {
                log_warning("too many errors to repair: interleave:%d count:%d", k, outer_error_locations.size());
                log_warning("attempting erasure decoding");
                rso.repair_interleaved(&message[0], &rso_ecc[0], k * rsi.message_len, rsi.message_len);
            }

            outer_error_locations.clear();
        }

        // TODO: sanity check on updated ecc
        std::copy_n(&rso_ecc[0], rso.interleaved_ecc_len, &ecc[0]);
        rsi.encode_blocks(&message[0], rso.message_len * interleave, &rsi_ecc[0]);
        rsi.encode_blocks(&rso_ecc[0], rso.ecc_len * interleave, &rsio_ecc[0]);

        return false;
    }
};
