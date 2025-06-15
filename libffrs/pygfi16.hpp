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


template<typename GFT, typename GF>
class ntt_naive {
public:
    template<typename T, typename U>
    inline void ntt(const GFT root, const T input[], size_t input_size, U output[], size_t output_size) const {
        auto& gf = GF::cast(this);
        for (size_t i = 0; i < output_size; ++i) {
            GFT acc = 0;
            for (size_t j = 0; j < input_size; ++j) {
                GFT exp = gf.pow(root, i * j);
                GFT prod = gf.mul(input[j], exp);
                acc = gf.add(acc, prod);
            }
            output[i] = acc;
        }
    }
};


template<typename GFT, typename GF>
class ntt_ct_iter {
public:
    template<typename T, typename U>
    inline void ntt(const GFT root, const T input[], size_t input_size, U output[], size_t output_size) const {
        copy_rbo(input, input_size, output, output_size);
        ntt_butterfly(root, input, input_size, output, output_size);
    }

    template<typename T, typename U>
    inline void ntt_rbo(const GFT root, const T input[], size_t input_size, U output[], size_t output_size) const {
        if (output_size < input_size)
            throw std::runtime_error("Output buffer smaller than input");

        std::copy_n(input, input_size, output);
        ntt_butterfly(root, input, input_size, output, output_size);
    }

    template<typename T, typename U>
    inline void copy_rbo(const T input[], size_t input_size, U output[], size_t output_size) const {
        output[0] = input[0];

        if (output_size <= 2)
            return;

        size_t l, rev = 0;
        for (size_t i = 1; i < input_size; ++i) {
            for (l = output_size >> 1; rev + l >= output_size; l >>= 1);
            rev = (rev & (l - 1)) + l;
            // rev = rbo16(i) >> (16 - nbits);
            output[rev] = input[i];
        }
    }

private:
    template<typename T, typename U>
    inline void ntt_butterfly(const GFT root, const T /*input*/[], size_t /*input_size*/, U output[], size_t output_size) const {
        auto& gf = GF::cast(this);
        for (size_t stride = 1, exp_f = output_size >> 1; stride < output_size; stride <<= 1, exp_f >>= 1) {
            for (size_t start = 0; start < output_size; start += stride * 2) {
                // For each pair of the CT butterfly operation.
                for (size_t i = start; i < start + stride; ++i) {
                    // j = i - start
                    GFT w = gf.pow(root, exp_f * (i - start));

                    // Cooley-Tukey butterfly
                    T a = output[i];
                    T b = output[i+stride];
                    GFT m = gf.mul(w, b);
                    output[i] = gf.add(a, m);
                    output[i+stride] = gf.sub(a, m);
                }
            }
        }
    }

    uint16_t rbo16(uint16_t b) const {
        b = ((b >> 1) & 0x5555) | ((b & 0x5555) << 1);
        b = ((b >> 2) & 0x3333) | ((b & 0x3333) << 2);
        b = ((b >> 4) & 0x0f0f) | ((b & 0x0f0f) << 4);
        b = ((b >> 8)         ) | ((b         ) << 8);
        return b;
    }
};


using GFi16 = ffrs::GF<uint32_t,
    ffrs::gf_data,
    // ffrs::gf_add_mod,
    // ffrs::gf_mul_mod,
    ffrs::gf_add_i16_shift,
    ffrs::gf_mul_i16_shift,
    ffrs::gf_exp_log_lut<ffrs::gf_mul_mod, 65537>::type,
    // ntt_naive,
    ntt_ct_iter
    >;


class PyGFi16 : public GFi16 {
public:
    inline PyGFi16(uint32_t prime, uint32_t primitive):
        GFi16(gf_data(prime, 1, primitive, 0))
    { }

    inline py::bytearray py_ntt(buffer_ro<GFT> buf, GFT root) {
        if (!ffrs::detail::is_power_of_two(buf.size))
            throw py::value_error("Buffer size must be a power of 2: " + std::to_string(buf.size));

        auto res =  py::bytearray(nullptr, buf.size * sizeof(uint32_t));

        auto data = PyByteArray_AsString(res.ptr());
        GFT *data_ptr = reinterpret_cast<GFT *>(data);

        ntt(root, buf.data, buf.size, data_ptr, buf.size);

        return res;
    }

    inline py::bytearray py_ntt16(buffer_ro<uint16_t> buf, GFT root) {
        if (!ffrs::detail::is_power_of_two(buf.size))
            throw py::value_error("Buffer size must be a power of 2: " + std::to_string(buf.size));

        auto res =  py::bytearray(nullptr, buf.size * sizeof(uint16_t));

        auto data = PyByteArray_AsString(res.ptr());
        uint16_t *data_ptr = reinterpret_cast<uint16_t *>(data);

        ntt(root, buf.data, buf.size, data_ptr, buf.size);

        return res;
    }

    inline py::bytearray py_ntt_blocks16(buffer_ro<uint16_t> buf, GFT root, size_t block_size) {
        if (!ffrs::detail::is_power_of_two(block_size))
            throw py::value_error("Block size must be a power of 2: " + std::to_string(block_size));
        if (buf.size % block_size != 0)
            throw py::value_error("Buffer size must be a multiple of block_size: " + std::to_string(buf.size));

        auto res =  py::bytearray(nullptr, buf.size * sizeof(uint16_t));
        auto data = PyByteArray_AsString(res.ptr());
        uint16_t *data_ptr = reinterpret_cast<uint16_t *>(data);

        for (size_t offset = 0; offset < buf.size; offset += block_size)
            ntt(root, &buf[offset], block_size, &data_ptr[offset], block_size);

        return res;
    }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyGFi16>(m, "GFi16")
            .def_property_readonly("prime", [](PyGFi16& self) { return self.prime; }, ":math:`p`")
            .def_property_readonly("power", [](PyGFi16& self) { return self.power; }, "Always 1")
            .def_property_readonly("primitive", [](PyGFi16& self) { return self.primitive; }, "Primitive value used to generate the field")
            .def_property_readonly("field_elements", [](PyGFi16& self) { return self.field_elements; }, ":math:`p`")
            .def(py::init<uint32_t, uint32_t>(), R"(
                Instantiate type for operations over :math:`GF(p)`

                Args:
                    prime : :math:`p` -- prime order of the field
                    primitive : :math:`a` -- primitive value used to generate the field
                )",
                "prime"_a, "primitive"_a
            )
            .def("mul", static_cast<GFT (PyGFi16::*)(GFT const&, GFT const&) const>(&PyGFi16::mul), R"(Multiplication: :math:`\text{lhs} \times \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("add", static_cast<GFT (PyGFi16::*)(GFT const&, GFT const&) const>(&PyGFi16::add), R"(Addition: :math:`\text{lhs} + \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("sub", static_cast<GFT (PyGFi16::*)(GFT const&, GFT const&) const>(&PyGFi16::sub), R"(Subtraction: :math:`\text{lhs} - \text{rhs}`)", "lhs"_a, "rhs"_a)
            .def("inv", &PyGFi16::inv, R"(Reciprocal: :math:`\frac{1}{\text{value}}`)", "value"_a)
            .def("div", &PyGFi16::div, R"(Division: :math:`\frac{\text{num}}{\text{den}}`)", "num"_a, "den"_a)
            .def("exp", &PyGFi16::exp, R"(Exponential function: :math:`a^{\text{value}}`)", "value"_a)
            .def("log", &PyGFi16::log, R"(Logarithm: :math:`\log_a (\text{value})`)", "value"_a)
            .def("pow", &PyGFi16::pow, R"(Power: :math:`\text{base}^\text{exponent}`)", "base"_a, "exponent"_a)
            .def("ntt", cast_args(&PyGFi16::py_ntt), R"(Number-theoretic transform (NTT) on a buffer)", "buf"_a, "root"_a)
            .def("ntt16", cast_args(&PyGFi16::py_ntt16), R"(Number-theoretic transform (NTT) on a buffer)", "buf"_a, "root"_a)
            .def("ntt_blocks16", cast_args(&PyGFi16::py_ntt_blocks16), R"(Number-theoretic transform (NTT) on 16-bit buffer blocks)", "buf"_a, "root"_a, "block_size"_a)
            .doc() = R"(
            Finite-field operations for prime fields <= 65537
            )";
    }
};
