/**************************************************************************
 * pyGFi32.hpp
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

#pragma once

#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "reed_solomon.hpp"
#include "util.hpp"

namespace py = pybind11;


template<typename GFT, typename GF>
class ntt_naive {
public:
    template<typename T, typename U>
    inline void ntt(const GFT root, const T input[], size_t input_size, U output[]) const {
        auto& gf = GF::cast(this);
        for (size_t i = 0; i < input_size; ++i) {
            GFT acc = 0;
            for (size_t j = 0; j < input_size; ++j) {
                GFT exp = gf.pow(root, i * j);
                // GFT exp = 1;
                // for (size_t k = 0; k < i*j; ++k)
                //     exp = (exp * root) % gf.prime;

                GFT prod = gf.mul(input[j], exp);
                // GFT prod = (input[j] * exp) % gf.prime;

                acc = gf.add(acc, prod);
                // acc = (acc + prod) % gf.prime;
            }
            output[i] = acc;
        }
    }
};


using GFi32 = ffrs::GF<uint32_t,
    ffrs::gf_data,
    ffrs::gf_add_mod,

    ffrs::gf_mul_mod,
    ffrs::gf_exp_log_lut<ffrs::gf_mul_mod, 65537>::type,

    ntt_naive
    >;


class PyGFi32 : public GFi32 {
public:
    inline PyGFi32(uint32_t prime, uint32_t primitive):
        GFi32(gf_data(prime, 1, primitive, 0))
    { }

    uint8_t rbo8(uint8_t b) {
        return ((b * 0x80200802ull) & 0x0884422110ull) * 0x0101010101ull >> 32;
    }

    uint32_t rbo32(uint32_t b) {
        b = ((b >> 1) & 0x55555555) | ((b & 0x55555555) << 1);
        b = ((b >> 2) & 0x33333333) | ((b & 0x33333333) << 2);
        b = ((b >> 4) & 0x0f0f0f0f) | ((b & 0x0f0f0f0f) << 4);
        b = ((b >> 8) & 0x00ff00ff) | ((b & 0x00ff00ff) << 8);
        b = (b >> 16) | (b << 16);
        return b;
    }

    inline py::bytearray py_ntt(buffer_ro<uint32_t> buf, GFT root) {
        if (buf.size == 0 || (buf.size & (buf.size - 1)) != 0)
            throw py::value_error("Buffer size must be a power of two: " + std::to_string(buf.size));

        auto res =  py::bytearray(nullptr, buf.size * sizeof(uint32_t));

        auto data = PyByteArray_AsString(res.ptr());
        GFT *data_ptr = reinterpret_cast<GFT *>(data);

        ntt(root, buf.data, buf.size, data_ptr);

        return res;
    }

    inline py::bytearray py_ntt16(buffer_ro<uint16_t> buf, GFT root) {
        if (buf.size == 0 || (buf.size & (buf.size - 1)) != 0)
            throw py::value_error("Buffer size must be a power of two: " + std::to_string(buf.size));

        auto res =  py::bytearray(nullptr, buf.size * sizeof(uint16_t));

        auto data = PyByteArray_AsString(res.ptr());
        uint16_t *data_ptr = reinterpret_cast<uint16_t *>(data);

        ntt(root, buf.data, buf.size, data_ptr);

        return res;
    }

    inline py::bytearray py_ntt_blocks16(buffer_ro<uint16_t> buf, GFT root, size_t block_size) {
        if (block_size == 0 || (block_size & (block_size - 1)) != 0)
            throw py::value_error("Block size must be a power of two: " + std::to_string(block_size));
        if (buf.size % block_size != 0)
            throw py::value_error("Buffer size must be a multiple of block_size: " + std::to_string(buf.size));

        auto res =  py::bytearray(nullptr, buf.size * sizeof(uint16_t));

        auto data = PyByteArray_AsString(res.ptr());
        uint16_t *data_ptr = reinterpret_cast<uint16_t *>(data);

        for (size_t offset = 0; offset < buf.size; offset += block_size)
            ntt(root, buf.data + offset, block_size, data_ptr);

        return res;
    }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyGFi32>(m, "GFi32")
            .def_property_readonly("prime", [](PyGFi32& self) { return self.prime; }, "Always 2")
            .def_property_readonly("power", [](PyGFi32& self) { return self.power; }, "Always 8")
            .def_property_readonly("primitive", [](PyGFi32& self) { return self.primitive; }, "Primitive value used to generate the field")
            .def_property_readonly("poly1", [](PyGFi32& self) { return self.poly1; }, "Masked irreducible polynomial, excluding MSb")
            .def_property_readonly("field_elements", [](PyGFi32& self) { return self.field_elements; }, "Always 256")
            .def(py::init<uint32_t, uint32_t>(), R"(
                Instantiate type for operations over :math:`GF(p^n)/P`

                Args:
                    prime : :math:`p` -- prime order of the field
                    primitive : :math:`a` -- primitive value used to generate the field
                )",
                "prime"_a, "primitive"_a
            )
            .def("mul", &PyGFi32::mul, R"(Multiplication: :math:`\text{lhs} \times \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("add", &PyGFi32::add, R"(Addition: :math:`\text{lhs} + \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("sub", &PyGFi32::sub, R"(Subtraction: :math:`\text{lhs} - \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("inv", &PyGFi32::inv, R"(Reciprocal: :math:`\frac{1}{\text{value}}`)", "value"_a)
            .def("div", &PyGFi32::div, R"(Division: :math:`\frac{\text{num}}{\text{den}}`)", "num"_a, "den"_a)
            .def("exp", &PyGFi32::exp, R"(Exponential function: :math:`a^{\text{value}}`)", "value"_a)
            .def("log", &PyGFi32::log, R"(Logarithm: :math:`\log_a (\text{value})`)", "value"_a)
            .def("pow", &PyGFi32::pow, R"(Power: :math:`\text{base}^\text{exponent}`)", "base"_a, "exponent"_a)
            .def("ntt", cast_args(&PyGFi32::py_ntt), R"(Number-theoretic transform (NTT) on a buffer)", "buf"_a, "root"_a)
            .def("ntt_blocks16", cast_args(&PyGFi32::py_ntt_blocks16), R"(Number-theoretic transform (NTT) on 16-bit buffer blocks)", "buf"_a, "root"_a, "block_size"_a)
            .doc() = R"(
                Finite-field operations for prime fields <= 65537
            )";
    }
};
