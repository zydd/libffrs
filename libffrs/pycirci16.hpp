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
    PyRSi16md rsoc;

public:
    const size_t message_len;
    const size_t rsi_chunk_ecc_len;
    const size_t rso_chunk_ecc_len;
    const size_t rsoc_chunk_ecc_len;
    const size_t ecc_len;

    inline PyCIRCi16(
        size_t inner_block_len,
        size_t inner_ecc_len,
        size_t outer_block_len,
        size_t outer_ecc_len,
        size_t outer_interleave,
        std::optional<bool> simd_x4,
        std::optional<bool> simd_x8,
        std::optional<bool> simd_x16
    ) :
        rsi(inner_block_len, {}, inner_ecc_len, 3, 1, simd_x4, simd_x8, simd_x16),
        rso(outer_block_len, {}, outer_ecc_len, 3, rsi.message_len * outer_interleave, simd_x4, simd_x8, simd_x16),
        rsoc(inner_block_len, {}, inner_ecc_len, 3, outer_ecc_len, simd_x4, simd_x8, simd_x16),
        message_len(rso.chunk_len),
        rsi_chunk_ecc_len(rsi.ecc_len * rso.message_len * outer_interleave),
        rso_chunk_ecc_len(rso.ecc_len * rso.interleave),
        rsoc_chunk_ecc_len(rsoc.ecc_len * rsoc.interleave),
        ecc_len(rsi_chunk_ecc_len + rso_chunk_ecc_len + rsoc_chunk_ecc_len)
    { }

    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        py_assert(buf.size % rso.chunk_len == 0, std::to_string(buf.size));
        py_assert(buf.size % rsi.message_len == 0, std::to_string(buf.size));
        size_t rsi_blocks = buf.size / rsi.message_len;
        py_assert(rsi_blocks * rsi.message_len == buf.size, std::to_string(buf.size));

        auto output = py::bytearray(nullptr, ecc_len * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));
        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{rsi.vec_align}) uint32_t[rso_chunk_ecc_len]);

        rsi.encode_full_blocks(&buf[0], rsi_blocks, 0, &output_data[0]);
        rso.encode_chunk(&buf[0], &temp[0]);
        std::copy_n(&temp[0], rso_chunk_ecc_len, &output_data[rsi_chunk_ecc_len]);
        rsoc.encode_chunk(&temp[0], &output_data[rsi_chunk_ecc_len + rso_chunk_ecc_len]);

        return output;
    }

    inline bool py_repair(buffer_rw<uint16_t> message, buffer_rw<uint16_t> ecc) {
        py_assert(message.size == message_len, std::to_string(message.size));
        py_assert(ecc.size == ecc_len, std::to_string(ecc.size));

        // TODO: handle outer_interleave
        size_t inner_blocks = message.size / rsi.message_len;
        py_assert(inner_blocks == rso.message_len);

        auto temp = std::unique_ptr<uint32_t[]>(new(std::align_val_t{rsi.vec_align}) uint32_t[rsi.block_len]);

        py_assert(rsi.ecc_len == rsoc.ecc_len);
        auto synds = std::unique_ptr<uint32_t[]>(new(std::align_val_t{rsi.vec_align}) uint32_t[rso.block_len * rsi.ecc_len]);

        rsi.synd_blocks(&message[0], &ecc[0], rso.message_len, &temp[0], &synds[0]);
        rsoc.synd_chunk(&ecc[rsi_chunk_ecc_len], &ecc[rsi_chunk_ecc_len + rso_chunk_ecc_len], &temp[0], &synds[rso.message_len * rsi.ecc_len]);

        std::vector<size_t> outer_error_locations;
        for (size_t i = 0; i < rso.block_len; ++i) {
            for (size_t j = 0; j < rsi.ecc_len; ++j) {
                if (synds[i * rsi.ecc_len + j] != 0) {
                    outer_error_locations.push_back(i);
                    break;
                }
            }
        }

        // py::print("\nerr locations:", outer_error_locations, "/", rso.block_len);

        // // TODO: repair errors in RSO ECC
        // for (size_t k = outer_error_locations.size() - 1; outer_error_locations.size() && outer_error_locations[k] >= rso.message_len; --k) {
        //     size_t j = rsi_chunk_ecc_len + outer_error_locations[k] - rso.message_len;
        //     std::vector<size_t> zeros;
        //     for (size_t i = 0; i < rso.interleave; ++i) {
        //         if (ecc[j + i * rso.ecc_len] == 0) {
        //             zeros.push_back(i);
        //         }
        //     }

        //     py::print("zeros:", zeros);
        //     py::print("synds:", std::vector<uint16_t>(&synds[outer_error_locations[k] * rsi.ecc_len], &synds[outer_error_locations[k] * rsi.ecc_len + rsi.ecc_len]));
        // }

        if (outer_error_locations.empty())
            return true;

        rso.repair_chunk(&message[0], &ecc[rsi_chunk_ecc_len], outer_error_locations);

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

            .def_property_readonly("message_len", [](PyCIRCi16& self) { return self.message_len; })
            .def_property_readonly("message_size", [](PyCIRCi16& self) { return self.message_len * sizeof(uint16_t); })
            .def_property_readonly("ecc_len", [](PyCIRCi16& self) { return self.ecc_len; })
            .def_property_readonly("ecc_size", [](PyCIRCi16& self) { return self.ecc_len * sizeof(uint16_t); })

            .def("encode", cast_args(&PyCIRCi16::py_encode), R"(Encode data)", "buffer"_a)
            .def("repair", cast_args(&PyCIRCi16::py_repair), R"(Repair data)", "message"_a, "ecc"_a);
    }
};
