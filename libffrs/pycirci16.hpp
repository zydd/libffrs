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

public:
    const size_t message_len;
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
        message_len(rso.chunk_len),
        ecc_len(rsi.ecc_len * rso.message_len + rso.ecc_len * rso.interleave)
    { }

    inline py::bytearray py_encode(buffer_ro<uint16_t> buf) {
        py_assert(buf.size % rso.chunk_len == 0, std::to_string(buf.size));
        py_assert(buf.size % rsi.message_len == 0, std::to_string(buf.size));
        size_t blocks = buf.size / rsi.message_len;

        auto output = py::bytearray(nullptr, ecc_len * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));

        rsi.encode_full_blocks(&buf[0], blocks, &output_data[0]);
        rsi.encode_chunk(&buf[0], &output_data[blocks * rsi.ecc_len]);

        return output;
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

            .def("encode", cast_args(&PyCIRCi16::py_encode), R"(Encode data)", "buffer"_a);
    }
};
