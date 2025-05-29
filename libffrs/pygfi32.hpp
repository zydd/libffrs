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


using GFi32 = ffrs::GF<uint32_t,
    ffrs::gf_data,
    ffrs::gf_add_mod,

    ffrs::gf_mul_mod,
    ffrs::gf_exp_log_lut<ffrs::gf_mul_mod, 65536>::type
    >;


class PyGFi32 : public GFi32 {
public:
    inline PyGFi32(uint32_t prime, uint32_t primitive):
        GFi32(gf_data(prime, 1, primitive, 0))
    { }


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
            .doc() = R"(
                Finite-field operations optimized for :math:`GF(2^8)`
            )";
    }
};
