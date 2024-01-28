/**************************************************************************
 * pyffrs.cpp
 *
 * Copyright 2024 Gabriel Machado
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

#include <pybind11/pybind11.h>

#include <libffrs/reed_solomon.hpp>

#include "util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;


using GF256 = ffrs::GF<uint8_t,
    ffrs::gf_add_xor,

    ffrs::gf_exp_log_lut<ffrs::gf_mul_cpu_pw2, 256>::type,

    // ffrs::gf_mul_cpu_pw2,
    // ffrs::gf_mul_lut<ffrs::gf_mul_cpu_pw2, 256>::type,
    ffrs::gf_mul_exp_log_lut,

    ffrs::gf_wide_mul<uint64_t>::type,
    ffrs::gf_poly_deriv_pw2,
    ffrs::gf_poly
    >;

template<typename GF>
using RS256 = ffrs::RS<GF,
    ffrs::rs_generator<256>::type,

    // ffrs::rs_encode_basic,
    // ffrs::rs_encode_basic_v2,
    // ffrs::rs_encode_lut_pw2<256>::type,
    // ffrs::rs_encode_slice_pw2<uint64_t, 8, 256, 16>::type,
    ffrs::rs_encode_slice_pw2_dispatch,

    // ffrs::rs_synds_basic<256>::type,
    ffrs::rs_synds_lut_pw2<uint32_t, 255>::type,

    // ffrs::rs_roots_eval_basic,
    // ffrs::rs_roots_eval_uint8_chien,
    ffrs::rs_roots_eval_lut_pw2<uint64_t>::type,

    ffrs::rs_decode
    >;


class PyGF256 : public GF256 {
public:
    inline PyGF256(GFT prime, GFT power, GFT primitive, GFT poly1):
        GF256(prime, power, primitive, poly1)
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
        if (buf2.size != sizeof(uint64_t))
            throw py::value_error("x must have 8 bytes");
        uint64_t res = poly_eval_wide(buf1.data, buf1.size, buf2.data[0]);
        return py::bytearray(reinterpret_cast<const char *>(&res), sizeof(uint64_t));
    }
    inline py::bytearray py_mul8(buffer_ro<uint64_t> buf1, buffer_ro<uint64_t> buf2) {
        if (buf1.size != sizeof(uint64_t) || buf2.size != sizeof(uint64_t))
            throw py::value_error("Can only multiply 8-byte wide arrays");
        uint64_t res = mul_wide(buf1.data[0], buf2.data[0]);
        return py::bytearray(reinterpret_cast<const char *>(&res), sizeof(uint64_t));
    }
};


class PyRS256 : public RS256<PyGF256> {
public:
    inline PyRS256(uint8_t ecc_len):
        RS256<PyGF256>(PyGF256(2, 8, 2, 0x1d), ecc_len)
    { }

    inline void py_encode(buffer_rw<uint8_t> buf) {
        size_t msg_size = buf.size - ecc_len;
        encode(buf.data, msg_size, &buf.data[msg_size]);
    }

    inline size_t compute_output_size(size_t input_size, size_t block_size) {
        size_t full_blocks = input_size / block_size;
        size_t remainder = input_size % block_size;

        size_t output_size = full_blocks * (block_size + this->ecc_len);

        if (remainder > 0)
            output_size += remainder + this->ecc_len;

        return output_size;
    }

    inline py::bytearray py_encode_blocks(buffer_ro<uint8_t> buf, size_t input_block_size) {
        // encode_blocks is 10x slower in MSVC compared to gcc/clang if we call rs_encode_slice_pw2_dispatch::encode
        // performing the dispatch 1 level higher seems to do the trick, MSVC is still slower, but comparable
        switch (ecc_len) {
        case 2: return _encode_blocks(&ffrs::rs_encode_slice_2<PyGF256, RS256<PyGF256>>::encode, std::move(buf), input_block_size);
        case 4: return _encode_blocks(&ffrs::rs_encode_slice_4<PyGF256, RS256<PyGF256>>::encode, std::move(buf), input_block_size);
        case 8: return _encode_blocks(&ffrs::rs_encode_slice_8<PyGF256, RS256<PyGF256>>::encode, std::move(buf), input_block_size);
        default: return _encode_blocks(&ffrs::rs_encode_basic_v2<PyGF256, RS256<PyGF256>>::encode, std::move(buf), input_block_size);
        }
    }

    inline bool py_decode(buffer_rw<uint8_t> buf) {
        size_t msg_size = buf.size - ecc_len;
        return decode(buf.data, msg_size, &buf.data[msg_size]);
    }

private:
    template<typename Method>
    inline py::bytearray _encode_blocks(Method method, buffer_ro<uint8_t> buf, size_t input_block_size) {
        if (buf.size == 0 || input_block_size == 0)
            return {};

        size_t full_blocks = buf.size / input_block_size;

        size_t output_size = full_blocks * (input_block_size + this->ecc_len);

        // Last block will be smaller if input size is not divisible by input_block_size
        size_t input_remainder = buf.size - full_blocks * input_block_size;
        if (input_remainder > 0)
            output_size += input_remainder + this->ecc_len;

        size_t output_block_size = input_block_size + this->ecc_len;

        auto output = py::bytearray();
        PyByteArray_Resize(output.ptr(), output_size);
        assert(size_t(PyByteArray_Size(output.ptr())) == output_size);
        auto output_data = reinterpret_cast<uint8_t *>(PyByteArray_AsString(output.ptr()));

        for (size_t block = 0; block < buf.size / input_block_size; ++block) {
            auto output_block = &output_data[block * output_block_size];
            std::copy_n(&buf.data[block * input_block_size], input_block_size, output_block);
            (this->*method)(output_block, input_block_size, &output_block[input_block_size]);
        }

        if (input_remainder > 0) {
            auto output_block = &output_data[output_size - input_remainder - this->ecc_len];
            std::copy_n(&buf.data[buf.size - input_remainder], input_remainder, output_block);
            (this->*method)(output_block, input_remainder, &output_block[input_remainder]);
        }

        return output;
    }

};


PYBIND11_MODULE(ffrs, m) {
    m.attr("__version__") = VERSION_INFO;
    m.doc() = R"(
        FFRS main module
    )";

    py::class_<PyGF256>(m, "GF256")
        .def_property_readonly("prime", [](PyGF256& self) { return self.prime; })
        .def_property_readonly("power", [](PyGF256& self) { return self.power; })
        .def_property_readonly("primitive", [](PyGF256& self) { return self.primitive; })
        .def_property_readonly("poly1", [](PyGF256& self) { return self.poly1; })
        .def_property_readonly("field_elements", [](PyGF256& self) { return self.field_elements; })
        .def(py::init<uint8_t, uint8_t, uint8_t, uint8_t>(), R"(
            Instantiate type for operations over :math:`GF(p^n)/P`

            Args:
                prime : :math:`p` -- always 2
                power : :math:`n` -- always 8
                primitive : :math:`a` -- primitive value used to generate the field
                polynomial : :math:`P` -- irreducible polynomial used to generate the field
            )",
            "prime"_a = 2, "power"_a = 8, "primitive"_a = 2, "poly1"_a = 0x1d
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

    py::class_<PyRS256>(m, "RS256")
        .def_property_readonly("ecc_len", [](PyRS256& self) { return self.ecc_len; })
        .def_property_readonly("gf", [](PyRS256& self) { return self.gf; })
        .def_property_readonly("generator", [](PyRS256& self) {
            return py::bytes(reinterpret_cast<const char *>(self.generator), self.ecc_len + 1); })
        .def(py::init<uint8_t>(), R"()", "ecc_len"_a)
        .def("__sizeof__", [](PyRS256& self) { return sizeof(self); })
        .def("encode", cast_args(&PyRS256::py_encode), R"(Systematic encode)", "buffer"_a)
        .def("encode_blocks", cast_args(&PyRS256::py_encode_blocks), R"(Encode blocks)", "buffer"_a, "input_block_size"_a)
        .def("decode", cast_args(&PyRS256::py_decode), R"(Systematic decode)", "buffer"_a)
        .doc() = R"(Reed-Solomon coding over :math:`GF(2^8)`)";
}
