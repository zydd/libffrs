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

# pragma once

#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
    const size_t outer_interleave;

    inline PyCIRC16(
        size_t inner_block_len,
        size_t inner_ecc_len,
        size_t outer_block_len,
        size_t outer_ecc_len,
        size_t outer_interleave,
        size_t primitive,
        std::optional<bool> simd_x4,
        std::optional<bool> simd_x8,
        std::optional<bool> simd_x16
    ):
        rsi(inner_block_len, inner_ecc_len, 1, primitive, simd_x4, simd_x8, simd_x16),
        rso(outer_block_len, outer_ecc_len, rsi.message_len * outer_interleave, primitive, simd_x4, simd_x8, simd_x16),
        block_len(inner_block_len * outer_block_len * outer_interleave),
        message_len(rso.interleaved_message_len),
        rsi_interleaved_ecc_len(rsi.ecc_len * rso.message_len * outer_interleave),
        rsio_ecc_len(rsi.ecc_len * rso.ecc_len * outer_interleave),
        ecc_len(rso.interleaved_ecc_len + rsi_interleaved_ecc_len + rsio_ecc_len),
        outer_interleave(outer_interleave)
    { }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyCIRC16>(m, "CIRC16")
            .def(py::init<
                    size_t,  // inner_block_len
                    size_t,  // inner_ecc_len
                    size_t,  // outer_block_len
                    size_t,  // outer_ecc_len
                    size_t,  // outer_interleave
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
                "outer_interleave"_a = 1,
                py::kw_only(),
                "primitive"_a = 3,
                "simd_x4"_a = py::none(),
                "simd_x8"_a = py::none(),
                "simd_x16"_a = py::none()
            )

            .def_property_readonly("rsi", [](PyCIRC16& self) -> auto const& { return self.rsi; })
            .def_property_readonly("rso", [](PyCIRC16& self) -> auto const& { return self.rso; })

            .def_property_readonly("block_len", [](PyCIRC16& self) { return self.block_len; })
            .def_property_readonly("block_size", [](PyCIRC16& self) { return self.block_len * sizeof(uint16_t); })
            .def_property_readonly("message_len", [](PyCIRC16& self) { return self.message_len; })
            .def_property_readonly("message_size", [](PyCIRC16& self) { return self.message_len * sizeof(uint16_t); })
            .def_property_readonly("ecc_len", [](PyCIRC16& self) { return self.ecc_len; })
            .def_property_readonly("ecc_size", [](PyCIRC16& self) { return self.ecc_len * sizeof(uint16_t); })
            .def_property_readonly("inner_block_len", [](PyCIRC16& self) { return self.rsi.block_len; })
            .def_property_readonly("inner_block_size", [](PyCIRC16& self) { return self.rsi.block_len * sizeof(uint16_t); })
            .def_property_readonly("inner_message_len", [](PyCIRC16& self) { return self.rsi.message_len; })
            .def_property_readonly("inner_message_size", [](PyCIRC16& self) { return self.rsi.message_len * sizeof(uint16_t); })
            .def_property_readonly("inner_ecc_len", [](PyCIRC16& self) { return self.rsi.ecc_len; })
            .def_property_readonly("inner_ecc_size", [](PyCIRC16& self) { return self.rsi.ecc_len * sizeof(uint16_t); })
            .def_property_readonly("outer_block_len", [](PyCIRC16& self) { return self.rso.block_len; })
            .def_property_readonly("outer_message_len", [](PyCIRC16& self) { return self.rso.message_len; })
            .def_property_readonly("outer_ecc_len", [](PyCIRC16& self) { return self.rso.ecc_len; })
            .def_property_readonly("outer_interleave", [](PyCIRC16& self) { return self.outer_interleave; })

            .def("encode", cast_args(&PyCIRC16::py_encode), R"(Encode data)", "buffer"_a)
            .def("repair", cast_args(&PyCIRC16::py_repair), R"(Repair data)", "message"_a, "ecc"_a);
    }

private:
    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        py_assert(buf.size % rso.interleaved_message_len == 0, std::to_string(buf.size));
        py_assert(buf.size % rsi.message_len == 0, std::to_string(buf.size));
        size_t rsi_blocks = buf.size / rsi.message_len;
        py_assert(rsi_blocks * rsi.message_len == buf.size, std::to_string(buf.size));

        auto output = py::bytearray(nullptr, ecc_len * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));
        auto temp = new_aligned<GFT>(rso.interleaved_ecc_len, rsi.vec_align);

        auto rso_ecc = &output_data[0];  // size = rso.interleaved_ecc_len
        auto rsi_ecc = &output_data[rso.interleaved_ecc_len];  // size = rsi_interleaved_ecc_len
        auto rsio_ecc = &output_data[rso.interleaved_ecc_len + rsi_interleaved_ecc_len];  // size = rsio_ecc_len

        rso.encode_interleaved(&buf[0], &temp[0]);
        std::copy_n(&temp[0], rso.interleaved_ecc_len, &rso_ecc[0]);
        rsi.encode_blocks(&buf[0], rsi_blocks, 0, &rsi_ecc[0]);
        rsi.encode_blocks(&temp[0], rso.ecc_len * outer_interleave, 0, &rsio_ecc[0]);

        return output;
    }

    inline bool py_repair(buffer_rw<uint16_t> message, buffer_rw<uint16_t> ecc) {
        py_assert(message.size == message_len, std::to_string(message.size));
        py_assert(ecc.size == ecc_len, std::to_string(ecc.size));

        size_t inner_blocks = message.size / rsi.message_len;
        py_assert(inner_blocks == rso.message_len * outer_interleave);

        auto buf_rsi = new_aligned<GFT>(rsi.block_len, rsi.vec_align);
        auto temp2_rsi = new_aligned<GFT>(rsi.block_len * 2, rsi.vec_align);
        auto synds_rsi = new_aligned<GFT>(rso.block_len * rsi.ecc_len * outer_interleave, rsi.vec_align);

        // Widen rso ecc to u32
        auto rso_ecc = new_aligned<GFT>(rso.interleaved_ecc_len, rsi.vec_align);
        std::copy_n(&ecc[0], rso.interleaved_ecc_len, &rso_ecc[0]);

        auto rsi_ecc = &ecc[rso.interleaved_ecc_len];
        auto rsio_ecc = &ecc[rso.interleaved_ecc_len + rsi_interleaved_ecc_len];

        rsi.synd_blocks(
            &message[0],
            &rsi_ecc[0],
            rso.message_len * outer_interleave,
            &temp2_rsi[0],
            &synds_rsi[0]
        );
        rsi.synd_blocks(
            &rso_ecc[0],
            &rsio_ecc[0],
            rso.ecc_len * outer_interleave,
            &temp2_rsi[0],
            &synds_rsi[rso.message_len * outer_interleave * rsi.ecc_len]
        );

        // Repair zeroes in rso ecc
        for (size_t k = 0; k < outer_interleave; ++k) {
            for (size_t i = 0; i < rso.ecc_len; ++i) {
                std::vector<size_t> inner_zero_locations;
                size_t rso_offset = (k + i * outer_interleave) * rsi.message_len;

                for (size_t j = 0; j < rsi.message_len; ++j) {
                    if (rso_ecc[rso_offset + j] == 0) {
                        inner_zero_locations.push_back(j);
                    }
                }

                size_t rsio_offset = (k + i * outer_interleave) * rsi.ecc_len;
                for (size_t j = 0; j < rsi.ecc_len; ++j) {
                    if (rsio_ecc[rsio_offset + j] == 0) {
                        inner_zero_locations.push_back(rsi.message_len + j);
                    }
                }
                if (inner_zero_locations.empty() || inner_zero_locations.size() >= rsi.ecc_len) {
                    continue;
                }

                // py::print("zeroes:", k, i, inner_zero_locations);

                std::copy_n(&rso_ecc[rso_offset], rsi.message_len, &buf_rsi[0]);
                std::copy_n(&rsio_ecc[rsio_offset], rsi.ecc_len, &buf_rsi[rsi.message_len]);
                rsi.repair_block(&buf_rsi[0], inner_zero_locations, &temp2_rsi[0]);
                std::copy_n(&buf_rsi[0], rsi.message_len, &rso_ecc[rso_offset]);

                if (std::all_of(inner_zero_locations.begin(), inner_zero_locations.end(),
                    [&](size_t zero_pos) { return (rso_ecc[rso_offset + zero_pos] & 0xffff) == 0; }))
                {
                    size_t synd_offset = (k + (rso.message_len + i) * outer_interleave) * rsi.ecc_len;
                    std::fill_n(&synds_rsi[synd_offset], rsi.ecc_len, 0);
                }
            }
        }

        std::vector<size_t> outer_error_locations;
        for (size_t k = 0; k < outer_interleave; ++k) {
            for (size_t i = 0; i < rso.block_len; ++i) {
                // size_t synd_offset = k * rsi.ecc_len + i * rsi.ecc_len * outer_interleave;
                size_t synd_offset = (k + i * outer_interleave) * rsi.ecc_len;
                if (std::any_of(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len], [](auto v) { return v != 0; })) {
                    // Repair zeroes in rsi ecc
                    if (std::any_of(&rsi_ecc[synd_offset], &rsi_ecc[synd_offset + rsi.ecc_len], [](auto v) { return v == 0; })) {
                        std::copy_n(&synds_rsi[synd_offset], rsi.ecc_len, &temp2_rsi[0]);
                        rsi.ecc_mix(&temp2_rsi[0]);

                        if (std::all_of(&temp2_rsi[0], &temp2_rsi[rsi.ecc_len], [](auto v) { return (v & 0xffff) == 0; })) {
                            continue;
                        }
                    }

                    outer_error_locations.push_back(i);
                }
            }

            if (outer_error_locations.empty())
                continue;

            // py::print("err locations:", k, outer_error_locations, "/", rso.block_len);
            // for (size_t i : outer_error_locations) {
            //     size_t synd_offset = (k + i * outer_interleave) * rsi.ecc_len;
            //     size_t msg_offset = (k + i * outer_interleave) * rsi.message_len;
            //     py::print("message:", k, i, std::vector(&message[msg_offset], &message[msg_offset + rsi.message_len]));
            //     py::print("ecc:", k, i, std::vector(&ecc[synd_offset], &ecc[synd_offset + rsi.ecc_len]));
            //     py::print("synds:", k, i, std::vector(&synds_rsi[synd_offset], &synds_rsi[synd_offset + rsi.ecc_len]));
            // }
            if (outer_error_locations.size() <= rso.ecc_len)
                rso.repair_interleaved(&message[0], &rso_ecc[0], k * rsi.message_len, rsi.message_len, outer_error_locations);
            else
                rso.repair_interleaved(&message[0], &rso_ecc[0], k * rsi.message_len, rsi.message_len);

            outer_error_locations.clear();
        }

        // TODO: sanity check on updated ecc
        std::copy_n(&rso_ecc[0], rso.interleaved_ecc_len, &ecc[0]);
        rsi.encode_blocks(&message[0], rso.message_len * outer_interleave, 0, &rsi_ecc[0]);
        rsi.encode_blocks(&rso_ecc[0], rso.ecc_len * outer_interleave, 0, &rsio_ecc[0]);

        return false;
    }
};
