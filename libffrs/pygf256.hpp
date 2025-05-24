/**************************************************************************
 * pygf256.hpp
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


using GF256 = ffrs::GF<uint8_t,
    ffrs::gf_data,
    ffrs::gf_add_xor,

    ffrs::gf_exp_log_lut<ffrs::gf_mul_cpu_pw2, 256>::type,

    // ffrs::gf_mul_cpu_pw2,
    // ffrs::gf_mul_lut<ffrs::gf_mul_cpu_pw2, 256>::type,
    ffrs::gf_mul_exp_log_lut,

    ffrs::gf_wide_mul<uint64_t>::type,
    ffrs::gf_poly_deriv_pw2,
    ffrs::gf_poly
    >;


class PyGF256 : public GF256 {
public:
    inline PyGF256(uint8_t primitive, uint16_t poly1):
        GF256(gf_data(2, 8, primitive, poly1 & 0xff))
    { }

    inline py::bytearray py_poly_add(buffer_ro<uint8_t> buf1, buffer_ro<uint8_t> buf2) {
        std::vector<uint8_t> output(std::max(buf1.size, buf2.size));
        poly_add(buf1.data, buf1.size, buf2.data, buf2.size, output.data());
        return py::bytearray(reinterpret_cast<const char *>(output.data()), output.size());
    }
    inline py::bytearray py_poly_sub(buffer_ro<uint8_t> buf1, buffer_ro<uint8_t> buf2) {
        std::vector<uint8_t> output(std::max(buf1.size, buf2.size));
        poly_sub(buf1.data, buf1.size, buf2.data, buf2.size, output.data());
        return py::bytearray(reinterpret_cast<const char *>(output.data()), output.size());
    }
    inline py::bytearray py_poly_mul(buffer_ro<uint8_t> buf1, buffer_ro<uint8_t> buf2) {
        std::vector<uint8_t> output(buf1.size + buf2.size - 1);
        poly_mul(buf1.data, buf1.size, buf2.data, buf2.size, output.data());
        return py::bytearray(reinterpret_cast<const char *>(output.data()), output.size());
    }
    inline py::bytearray py_poly_mod(buffer_ro<uint8_t> buf1, buffer_ro<uint8_t> buf2) {
        std::vector<uint8_t> output(buf2.size - 1);
        auto rem_size = poly_mod(buf1.data, buf1.size, buf2.data, buf2.size, output.data());
        return py::bytearray(reinterpret_cast<const char *>(output.data()), rem_size);
    }
    inline py::tuple py_poly_divmod(buffer_ro<uint8_t> buf1, buffer_ro<uint8_t> buf2) {
        std::vector<uint8_t> output(buf1.data, buf1.data + buf1.size);
        auto sep = ex_synth_div(output.data(), output.size(), buf2.data, buf2.size);
        return py::make_tuple(
            py::bytearray(reinterpret_cast<const char *>(&output[0]), sep),
            py::bytearray(reinterpret_cast<const char *>(&output[sep]), output.size() - sep));
    }
    inline py::bytearray py_poly_mod_x_n(buffer_ro<uint8_t> buf1, buffer_ro<uint8_t> buf2) {
        std::vector<uint8_t> output(buf2.size);
        poly_mod_x_n(buf1.data, buf1.size, buf2.data, buf2.size, output.data());
        return py::bytearray(reinterpret_cast<const char *>(output.data()), buf2.size);
    }
    inline uint8_t py_poly_eval(buffer_ro<uint8_t> buf, uint8_t x) {
        return poly_eval(buf.data, buf.size, x);
    }
    inline py::bytearray py_poly_eval8(buffer_ro<uint8_t> buf1, buffer_ro<uint64_t> buf2) {
        if (buf2.size > sizeof(uint64_t))
            throw py::value_error("can evaluate polynomial 8 points at most");
        uint64_t res = poly_eval_wide(buf1.data, buf1.size, buf2.data[0]);
        return py::bytearray(reinterpret_cast<const char *>(&res), std::min(buf2.size, sizeof(uint64_t)));
    }
    inline py::bytearray py_mul8(buffer_ro<uint64_t> buf1, buffer_ro<uint64_t> buf2) {
        if (buf1.size != sizeof(uint64_t) || buf2.size != sizeof(uint64_t))
            throw py::value_error("Can only multiply 8-byte wide arrays");
        uint64_t res = mul_wide(buf1.data[0], buf2.data[0]);
        return py::bytearray(reinterpret_cast<const char *>(&res), sizeof(uint64_t));
    }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyGF256>(m, "GF256")
            .def_property_readonly("prime", [](PyGF256& self) { return self.prime; }, "Always 2")
            .def_property_readonly("power", [](PyGF256& self) { return self.power; }, "Always 8")
            .def_property_readonly("primitive", [](PyGF256& self) { return self.primitive; }, "Primitive value used to generate the field")
            .def_property_readonly("poly1", [](PyGF256& self) { return self.poly1; }, "Masked irreducible polynomial, excluding MSb")
            .def_property_readonly("field_elements", [](PyGF256& self) { return self.field_elements; }, "Always 256")
            .def(py::init<uint8_t, uint16_t>(), R"(
                Instantiate type for operations over :math:`GF(p^n)/P`

                Args:
                    primitive : :math:`a` -- primitive value used to generate the field
                    polynomial : :math:`P` -- irreducible polynomial used to generate the field
                )",
                "primitive"_a = 2, "poly1"_a = 0x11d
            )
            .def("mul", &PyGF256::mul, R"(Multiplication: :math:`\text{lhs} \times \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("add", &PyGF256::add, R"(Addition: :math:`\text{lhs} + \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("sub", &PyGF256::sub, R"(Subtraction: :math:`\text{lhs} - \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("inv", &PyGF256::inv, R"(Reciprocal: :math:`\frac{1}{\text{value}}`)", "value"_a)
            .def("div", &PyGF256::div, R"(Division: :math:`\frac{\text{num}}{\text{den}}`)", "num"_a, "den"_a)
            .def("exp", &PyGF256::exp, R"(Exponential function: :math:`a^{\text{value}}`)", "value"_a)
            .def("log", &PyGF256::log, R"(Logarithm: :math:`\log_a (\text{value})`)", "value"_a)
            .def("pow", &PyGF256::pow, R"(Power: :math:`\text{base}^\text{exponent}`)", "base"_a, "exponent"_a)
            .def("poly_add", cast_args(&PyGF256::py_poly_add), R"(Add polynomials)", "p1"_a, "p2"_a)
            .def("poly_sub", cast_args(&PyGF256::py_poly_sub), R"(Subtract polynomials)", "p1"_a, "p2"_a)
            .def("poly_mul", cast_args(&PyGF256::py_poly_mul), R"(Multiply polynomials)", "p1"_a, "p2"_a)
            .def("poly_mod", cast_args(&PyGF256::py_poly_mod), R"(Polynomial remainder)", "p1"_a, "p2"_a)
            .def("poly_eval", cast_args(&PyGF256::py_poly_eval), R"(Evaluate polynomial at ``x``)", "poly"_a, "x"_a)
            .def("poly_eval8", cast_args(&PyGF256::py_poly_eval8),
                R"(Evaluate polynomial at 8 points on a single operation)", "poly"_a, "xs"_a)
            .def("mul8", cast_args(&PyGF256::py_mul8), R"(Multiply 8 values in a single operation)", "a"_a, "b"_a)
            .def("poly_divmod", cast_args(&PyGF256::py_poly_divmod), R"(
                    Polynomial quotient and remainder

                    Returns:
                        (quotient, remainder)
                )", "p1"_a, "p2"_a)
            .def("poly_mod_x_n", cast_args(&PyGF256::py_poly_mod_x_n), R"(
                    Shifted polynomial remainder

                    :math:`P \times X^n \mod (X^n + p_2)` where ``n = len(p2)``
                )", "p1"_a, "p2"_a)
            .doc() = R"(
                Finite-field operations optimized for :math:`GF(2^8)`
            )";
    }
};