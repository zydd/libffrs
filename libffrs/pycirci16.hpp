/**************************************************************************
 * pycirci16.hpp
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

#include "pyrsi16md.hpp"


namespace py = pybind11;


class PyCIRCi16 {
private:
    PyRSi16md rsi;
    PyRSi16md rso;
    // PyRSi16md rsoc;

public:
    const size_t block_len;
    const size_t message_len;
    const size_t rsi_chunk_ecc_len;
    const size_t rsi_rso_ecc_len;
    const size_t ecc_len;
    const size_t outer_interleave;

    inline PyCIRCi16(
        size_t inner_block_len,
        size_t inner_ecc_len,
        size_t outer_block_len,
        size_t outer_ecc_len,
        size_t outer_interleave,
        std::optional<bool> simd_x4,
        std::optional<bool> simd_x8,
        std::optional<bool> simd_x16
    ):
        rsi(inner_block_len, {}, inner_ecc_len, 3, 1, simd_x4, simd_x8, simd_x16),
        rso(outer_block_len, {}, outer_ecc_len, 3, rsi.message_len * outer_interleave, simd_x4, simd_x8, simd_x16),
        // rsoc(inner_block_len, {}, inner_ecc_len, 3, outer_ecc_len * outer_interleave, simd_x4, simd_x8, simd_x16),
        block_len(inner_block_len * outer_block_len * outer_interleave),
        message_len(rso.chunk_len),
        rsi_chunk_ecc_len(rsi.ecc_len * rso.message_len * outer_interleave),
        rsi_rso_ecc_len(rsi.ecc_len * rso.ecc_len * outer_interleave),
        ecc_len(rsi_chunk_ecc_len + rso.chunk_ecc_len + rsi_rso_ecc_len),
        // ecc_len(rsi_chunk_ecc_len + rso.chunk_ecc_len + rsoc.chunk_ecc_len),
        outer_interleave(outer_interleave)
    { }

    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        py_assert(buf.size % rso.chunk_len == 0, std::to_string(buf.size));
        py_assert(buf.size % rsi.message_len == 0, std::to_string(buf.size));
        size_t rsi_blocks = buf.size / rsi.message_len;
        py_assert(rsi_blocks * rsi.message_len == buf.size, std::to_string(buf.size));

        auto output = py::bytearray(nullptr, ecc_len * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));
        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{rsi.vec_align}) uint32_t[rso.chunk_ecc_len]);

        rsi.encode_full_blocks(&buf[0], rsi_blocks, 0, &output_data[0]);
        rso.encode_chunk(&buf[0], &temp[0]);
        std::copy_n(&temp[0], rso.chunk_ecc_len, &output_data[rsi_chunk_ecc_len]);

        // py_assert(rso.chunk_ecc_len == rsoc.chunk_message_len);
        // rsoc.encode_chunk(&temp[0], &output_data[rsi_chunk_ecc_len + rso.chunk_ecc_len]);
        py_assert(rso.chunk_ecc_len == rsi.message_len * rso.ecc_len * outer_interleave);
        rsi.encode_full_blocks(&temp[0], rso.ecc_len * outer_interleave, 0, &output_data[rsi_chunk_ecc_len + rso.chunk_ecc_len]);

        return output;
    }

    inline bool py_repair(buffer_rw<uint16_t> message, buffer_rw<uint16_t> ecc) {
        py_assert(message.size == message_len, std::to_string(message.size));
        py_assert(ecc.size == ecc_len, std::to_string(ecc.size));

        // TODO: handle outer_interleave
        size_t inner_blocks = message.size / rsi.message_len;
        py_assert(inner_blocks == rso.message_len * outer_interleave);

        auto buf_rsi = std::unique_ptr<uint32_t[]>(new(std::align_val_t{rsi.vec_align}) uint32_t[rsi.block_len]);
        auto temp_rsi = std::unique_ptr<uint32_t[]>(new(std::align_val_t{rsi.vec_align}) uint32_t[rsi.block_len]);
        auto synds_rsi = std::unique_ptr<uint32_t[]>(new(std::align_val_t{rsi.vec_align}) uint32_t[rso.block_len * rsi.ecc_len * outer_interleave]);

        rsi.synd_blocks(&message[0], &ecc[0], rso.message_len * outer_interleave, &temp_rsi[0], &synds_rsi[0]);

        // py_assert(rsi.ecc_len == rsoc.ecc_len);
        // rsoc.synd_chunk(&ecc[rsi_chunk_ecc_len], &ecc[rsi_chunk_ecc_len + rso.chunk_ecc_len], &temp_rsi[0], &synds_rsi[rso.message_len * outer_interleave * rsi.ecc_len]);

        auto ecc_rso = std::unique_ptr<uint32_t[]>(new(std::align_val_t{rsi.vec_align}) uint32_t[rso.chunk_ecc_len]);

        std::copy_n(&ecc[rsi_chunk_ecc_len], rso.chunk_ecc_len, &ecc_rso[0]);
        rsi.synd_blocks(&ecc_rso[0], &ecc[rsi_chunk_ecc_len + rso.chunk_ecc_len], rso.ecc_len * outer_interleave, &temp_rsi[0], &synds_rsi[rso.message_len * outer_interleave * rsi.ecc_len]);

        // Repair zeroes in rso ecc
        for (size_t k = 0; k < outer_interleave; ++k) {
            for (size_t i = 0; i < rso.ecc_len; ++i) {
                std::vector<size_t> inner_zero_locations;
                size_t offset = (k + i * outer_interleave) * rsi.message_len;

                for (size_t j = 0; j < rsi.message_len; ++j) {
                    if (ecc_rso[offset + j] == 0) {
                        inner_zero_locations.push_back(j);
                    }
                }

                size_t ecc_offset = rsi_chunk_ecc_len + rso.chunk_ecc_len + (k + i * outer_interleave) * rsi.ecc_len;
                for (size_t j = 0; j < rsi.ecc_len; ++j) {
                    if (ecc[ecc_offset + j] == 0) {
                        inner_zero_locations.push_back(rsi.message_len + j);
                    }
                }
                if (inner_zero_locations.empty()) {
                    continue;
                }

                // py::print("zeroes:", k, i, inner_zero_locations);

                std::copy_n(&ecc_rso[offset], rsi.message_len, &buf_rsi[0]);
                std::copy_n(&ecc[ecc_offset], rsi.ecc_len, &buf_rsi[rsi.message_len]);
                rsi.repair_block(&buf_rsi[0], inner_zero_locations, &temp_rsi[0]);
                std::copy_n(&buf_rsi[0], rsi.message_len, &ecc_rso[offset]);

                for (size_t zero_pos : inner_zero_locations) {
                    if (zero_pos < rsi.message_len)
                        py_assert((ecc_rso[offset + zero_pos] & 0xffff) == 0, std::to_string(ecc_rso[offset + zero_pos]));
                }

                size_t synd_offset = (k + (rso.message_len + i) * outer_interleave) * rsi.ecc_len;
                std::fill_n(&synds_rsi[synd_offset], rsi.ecc_len, 0);
            }
        }

        std::vector<size_t> outer_error_locations;
        for (size_t k = 0; k < outer_interleave; ++k) {
            for (size_t i = 0; i < rso.block_len; ++i) {
                // size_t synd_offset = k * rsi.ecc_len + i * rsi.ecc_len * outer_interleave;
                size_t synd_offset = (k + i * outer_interleave) * rsi.ecc_len;
                if (std::any_of(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len], [](auto v) { return v != 0; })) {
                    // Repair zeroes in rsi ecc
                    if (std::any_of(&ecc[synd_offset], &ecc[synd_offset + rsi.ecc_len], [](auto v) { return v == 0; })) {
                        std::copy_n(&synds_rsi[synd_offset], rsi.ecc_len, &temp_rsi[0]);
                        rsi.ecc_mix(&temp_rsi[0]);

                        if (std::all_of(&temp_rsi[0], &temp_rsi[rsi.ecc_len], [](auto v) { return (v & 0xffff) == 0; })) {
                            continue;
                        }
                    }

                    outer_error_locations.push_back(i);
                }
            }

            if (outer_error_locations.empty())
                continue;

            py::print("err locations:", k, outer_error_locations, "/", rso.block_len);
            // for (size_t i : outer_error_locations) {
            //     size_t synd_offset = (k + i * outer_interleave) * rsi.ecc_len;
            //     size_t msg_offset = (k + i * outer_interleave) * rsi.message_len;
            //     py::print("message:", k, i, std::vector(&message[msg_offset], &message[msg_offset + rsi.message_len]));
            //     py::print("ecc:", k, i, std::vector(&ecc[synd_offset], &ecc[synd_offset + rsi.ecc_len]));
            //     py::print("synds:", k, i, std::vector(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len]));
            // }
            rso.repair_chunk(&message[0], &ecc_rso[0], k * rsi.message_len, rsi.message_len, outer_error_locations);

            outer_error_locations.clear();
        }

        return false;
    }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyCIRCi16>(m, "CIRCi16")
            .def(py::init<
                    size_t,  // inner_block_len
                    size_t,  // inner_ecc_len
                    size_t,  // outer_block_len
                    size_t,  // outer_ecc_len
                    size_t,  // outer_interleave
                    std::optional<bool>,  // simd_x4
                    std::optional<bool>,  // simd_x8
                    std::optional<bool>   // simd_x16
                >(),
                R"(Cross-interleaved Reed-Solomon coder)",
                "inner_block_len"_a,
                "inner_ecc_len"_a,
                "outer_block_len"_a,
                "outer_ecc_len"_a,
                "outer_interleave"_a = 1,
                "simd_x4"_a = py::none(),
                "simd_x8"_a = py::none(),
                "simd_x16"_a = py::none()
            )

            .def_property_readonly("rsi", [](PyCIRCi16& self) -> auto const& { return self.rsi; })
            .def_property_readonly("rso", [](PyCIRCi16& self) -> auto const& { return self.rso; })

            .def_property_readonly("block_len", [](PyCIRCi16& self) { return self.block_len; })
            .def_property_readonly("block_size", [](PyCIRCi16& self) { return self.block_len * sizeof(uint16_t); })
            .def_property_readonly("message_len", [](PyCIRCi16& self) { return self.message_len; })
            .def_property_readonly("message_size", [](PyCIRCi16& self) { return self.message_len * sizeof(uint16_t); })
            .def_property_readonly("ecc_len", [](PyCIRCi16& self) { return self.ecc_len; })
            .def_property_readonly("ecc_size", [](PyCIRCi16& self) { return self.ecc_len * sizeof(uint16_t); })
            .def_property_readonly("outer_interleave", [](PyCIRCi16& self) { return self.outer_interleave; })

            .def("encode", cast_args(&PyCIRCi16::py_encode), R"(Encode data)", "buffer"_a)
            .def("repair", cast_args(&PyCIRCi16::py_repair), R"(Repair data)", "message"_a, "ecc"_a);
    }
};
