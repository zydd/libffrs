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

#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
    ffrs::rs_encode_slice_pw2_dispatch<255>::type,

    // ffrs::rs_synds_basic<256>::type,
    ffrs::rs_synds_lut_pw2<uint32_t, 255>::type,

    // ffrs::rs_roots_eval_basic,
    // ffrs::rs_roots_eval_uint8_chien,
    ffrs::rs_roots_eval_lut_pw2<uint64_t>::type,

    ffrs::rs_decode
    >;


class PyGF256 : public GF256 {
public:
    inline PyGF256(uint8_t primitive, uint16_t poly1):
        GF256(2, 8, primitive, poly1 & 0xff)
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
    size_t default_block_size = 255;

    inline PyRS256(
            std::optional<uint8_t> block_size,
            std::optional<uint8_t> message_len,
            std::optional<uint8_t> ecc_len,
            uint8_t primitive,
            uint16_t polynomial):
        PyRS256(_get_ecc_len(block_size, message_len, ecc_len),
                block_size.value_or(255),
                primitive, polynomial)
    {
        if (ecc_len && message_len) {
            if (block_size && *message_len + *ecc_len != *block_size) {
                throw py::value_error("block_len must be equal to message_len + ecc_len");
            } else if (!block_size) {
                set_default_block_size(size_t(*message_len) + size_t(*ecc_len));
            }
        }
    }

    inline PyRS256(uint8_t ecc_len, size_t block_size, uint8_t primitive, uint16_t polynomial):
        RS256<PyGF256>(PyGF256(primitive, polynomial), ecc_len)
    {
        set_default_block_size(block_size);
    }

    inline void set_default_block_size(size_t block_size) {
        if (block_size <= ecc_len)
            throw py::value_error("block_len must be greater than ecc_len");

        if (block_size > 255)
            throw py::value_error("block_len must be <= 255");

        default_block_size = block_size;
    }

    inline void py_encode(buffer_rw<uint8_t> buf) {
        size_t msg_size = buf.size - ecc_len;
        encode(buf.data, msg_size, &buf.data[msg_size]);
    }

    inline py::bytearray py_encode_blocks(buffer_ro<uint8_t> buf, std::optional<size_t> block_size) {
        size_t output_block_size = block_size.value_or(default_block_size);
        size_t input_block_size = output_block_size - ecc_len;

        if (buf.size == 0 || output_block_size == 0 || output_block_size <= ecc_len)
            return {};

        size_t full_blocks = buf.size / input_block_size;

        size_t output_size = full_blocks * (input_block_size + this->ecc_len);

        // Last block will be smaller if input size is not divisible by input_block_size
        size_t input_remainder = buf.size - full_blocks * input_block_size;
        if (input_remainder > 0)
            output_size += input_remainder + this->ecc_len;

        auto output = py::bytearray();
        PyByteArray_Resize(output.ptr(), output_size);
        assert(size_t(PyByteArray_Size(output.ptr())) == output_size);
        auto output_data = reinterpret_cast<uint8_t *>(PyByteArray_AsString(output.ptr()));

        for (size_t block = 0; block < buf.size / input_block_size; ++block) {
            auto output_block = &output_data[block * output_block_size];
            std::copy_n(&buf.data[block * input_block_size], input_block_size, output_block);
            encode(output_block, input_block_size, &output_block[input_block_size]);
        }

        if (input_remainder > 0) {
            auto output_block = &output_data[output_size - input_remainder - this->ecc_len];
            std::copy_n(&buf.data[buf.size - input_remainder], input_remainder, output_block);
            encode(output_block, input_remainder, &output_block[input_remainder]);
        }

        return output;
    }

    inline bool py_decode(buffer_rw<uint8_t> buf) {
        size_t msg_size = buf.size - ecc_len;
        return decode(buf.data, msg_size, &buf.data[msg_size]);
    }

    inline py::bytearray py_synds(buffer_ro<uint8_t> buf) {
        if (buf.size < ecc_len)
            return {};

        size_t msg_size = buf.size - ecc_len;
        synds_array_t synds_arr;

        synds(buf.data, msg_size, &buf.data[msg_size], synds_arr);

        return py::bytearray(reinterpret_cast<const char *>(synds_arr), ecc_len);
    }

private:
    inline uint8_t _get_ecc_len(
            std::optional<uint8_t> block_len,
            std::optional<uint8_t> message_len,
            std::optional<uint8_t> ecc_len) {
        if (ecc_len) {
            return *ecc_len;
        } else if (message_len && block_len) {
            if (*message_len >= *block_len) {
                throw py::value_error("block_len must be greater than message_len");
            }
            return block_len.value_or(default_block_size) - *message_len;
        } else {
            throw py::value_error("Must specify either (block_len, message_len) or ecc_len");
        }
    }
};


PYBIND11_MODULE(ffrs, m) {
    m.attr("__version__") = VERSION_INFO;
#if defined(__clang__)
    m.attr("compiler_info") = __VERSION__;
#elif defined(__GNUC__)
    m.attr("compiler_info") = "gcc" __VERSION__;
#elif defined(_MSC_VER)
    #define _str(a) _stringify(a)
    #define _stringify(a) #a
    m.attr("compiler_info") = "MSVC" _str(_MSC_VER);
#else
    m.attr("compiler_info") = "unknown";
#endif

    m.doc() = R"(
        FFRS main module
    )";

    py::class_<PyGF256>(m, "GF256")
        .def_property_readonly("prime", [](PyGF256& self) { return self.prime; })
        .def_property_readonly("power", [](PyGF256& self) { return self.power; })
        .def_property_readonly("primitive", [](PyGF256& self) { return self.primitive; })
        .def_property_readonly("poly1", [](PyGF256& self) { return self.poly1; })
        .def_property_readonly("field_elements", [](PyGF256& self) { return self.field_elements; })
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

    py::class_<PyRS256>(m, "RS256")
        .def_property_readonly("ecc_len", [](PyRS256& self) { return self.ecc_len; })

        .def_property("default_block_len",
            [](PyRS256& self) { return self.default_block_size; },
            &PyRS256::set_default_block_size)

        .def_property_readonly("default_block_msg_len",
            [](PyRS256& self) { return self.default_block_size - self.ecc_len; })

        .def_property_readonly("gf", [](PyRS256& self) -> auto const& { return self.gf; })

        .def_property_readonly("generator", [](PyRS256& self) {
            return py::bytes(reinterpret_cast<const char *>(self.generator), self.ecc_len + 1); })

        .def_property_readonly("generator_roots", [](PyRS256& self) {
            return py::bytes(reinterpret_cast<const char *>(self.generator_roots), self.ecc_len); })

        .def(py::init<std::optional<uint8_t>, std::optional<uint8_t>, std::optional<uint8_t>, uint8_t, uint16_t>(), R"(
            Instantiate a Reed-Solomon encoder with the given configuration

            Example:
                ``RS256(255, 223)``
                    Equivalent to ``RS256(ecc_len=32)``

                    Creates an encoder for 32 bytes of parity, capable of correcting up to 16 errors in a 255-byte block

            Args:
                block_len
                    | :math:`n` -- default block size used by :py:meth:`encode_blocks`
                    | ``block_len = message_len + ecc_len``

                message_len
                    | :math:`k` -- number of actual data bytes in a block
                    | Can be omitted if ``ecc_len`` is supplied

                ecc_len
                    | :math:`(n - k)` -- number of parity bytes in a block
                    |  Can be omitted if ``message_len`` is supplied

                primitive
                    :math:`a` -- primitive value for :py:class:`GF256`

                polynomial
                    :math:`P` -- irreducible polynomial for :py:class:`GF256`
            )",
            "block_len"_a = py::none(), "message_len"_a = py::none(), "ecc_len"_a = py::none(),
            "primitive"_a = 2, "polynomial"_a = 0x11d)

        .def("__sizeof__", [](PyRS256& self) { return sizeof(self); })

        .def("encode", cast_args(&PyRS256::py_encode),
            R"(Systematic encode)",
            "buffer"_a)

        .def("encode_blocks", cast_args(&PyRS256::py_encode_blocks), R"(
            Encode blocks

            .. note::
                The size of ``buffer`` should be a multiple of :py:attr:`RS256.default_block_msg_len`
                to allow concatenating the results of multiple calls to :py:meth:`encode_blocks`.
            )",
            "buffer"_a, "block_len"_a = py::none())

        .def("decode", cast_args(&PyRS256::py_decode),
            R"(Systematic decode)",
            "buffer"_a)

        .def("_synds", cast_args(&PyRS256::py_synds),
            R"(Compute syndromes)",
            "buffer"_a)

        .doc() = R"(Reed-Solomon coding over :math:`GF(2^8)`)";
}
