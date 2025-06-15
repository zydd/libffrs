/**************************************************************************
 * pyGFi16.hpp
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

#include "reed_solomon.hpp"
#include "util.hpp"

namespace py = pybind11;

using GFi16 = ffrs::GF<uint32_t,
    ffrs::gf_data,
    // ffrs::gf_add_mod,
    // ffrs::gf_mul_mod,
    ffrs::gf_add_i16_shift,
    ffrs::gf_mul_i16_shift,
    ffrs::gf_exp_log_lut<ffrs::gf_mul_mod, 0x10001>::type
    >;


class PyGFi16 : public GFi16 {
public:
    inline PyGFi16(uint32_t primitive):
        GFi16(gf_data(0x10001, 1, primitive, 0))
    { }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyGFi16>(m, "GFi16")
            .def_property_readonly("prime", [](PyGFi16& self) { return self.prime; }, ":math:`p`")
            .def_property_readonly("power", [](PyGFi16& self) { return self.power; }, "Always 1")
            .def_property_readonly("primitive", [](PyGFi16& self) { return self.primitive; }, "Primitive value used to generate the field")
            .def_property_readonly("field_elements", [](PyGFi16& self) { return self.field_elements; }, ":math:`p`")
            .def(py::init<uint32_t>(), R"(
                Instantiate type for operations over :math:`GF(65537)`

                Args:
                    primitive : :math:`a` -- primitive value used to generate the field
                )",
                "primitive"_a
            )
            .def("mul", static_cast<GFT (PyGFi16::*)(GFT const&, GFT const&) const>(&PyGFi16::mul), R"(Multiplication: :math:`\text{lhs} \times \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("add", static_cast<GFT (PyGFi16::*)(GFT const&, GFT const&) const>(&PyGFi16::add), R"(Addition: :math:`\text{lhs} + \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("sub", static_cast<GFT (PyGFi16::*)(GFT const&, GFT const&) const>(&PyGFi16::sub), R"(Subtraction: :math:`\text{lhs} - \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("inv", &PyGFi16::inv, R"(Reciprocal: :math:`\frac{1}{\text{value}}`)", "value"_a)
            .def("div", &PyGFi16::div, R"(Division: :math:`\frac{\text{num}}{\text{den}}`)", "num"_a, "den"_a)
            .def("exp", &PyGFi16::exp, R"(Exponential function: :math:`a^{\text{value}}`)", "value"_a)
            .def("log", &PyGFi16::log, R"(Logarithm: :math:`\log_a (\text{value})`)", "value"_a)
            .def("pow", &PyGFi16::pow, R"(Power: :math:`\text{base}^\text{exponent}`)", "base"_a, "exponent"_a)
            .doc() = R"(
            Finite-field operations for prime fields <= 65537
            )";
    }
};
