/**************************************************************************
 * pyntt.hpp
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

#include "ntt.hpp"


/**
 * Register NTT class into a pybind11 module
 * Cannot be included in the same translation unit as the instantiation of the templates in NTT
 */
inline void NTT::register_class(py::module &m) {
    using namespace pybind11::literals;

    py::class_<NTT>(m, "NTT")
        .def(py::init<PyGFi16 const&, size_t, size_t>(), R"(
            Args:
                gf: :class:`GFi16` instance
                block_len: Number of data symbols per block.
                ecc_len: Number of ECC symbols per block.
        )", "gf"_a, "block_len"_a, "ecc_len"_a)

        .def_property_readonly("ntt_len", [](NTT& self) { return self.ntt_len; })
        .def_property_readonly("ntt_size", [](NTT& self) { return self.ntt_len * sizeof(uint16_t); })
        .def_property_readonly("ntt4_size", [](NTT& self) { return self.ntt_len * 4 * sizeof(uint16_t); })
        .def_property_readonly("ntt8_size", [](NTT& self) { return self.ntt_len * 8 * sizeof(uint16_t); })
        .def_property_readonly("ntt16_size", [](NTT& self) { return self.ntt_len * 16 * sizeof(uint16_t); })
        .def_property_readonly("block_len", [](NTT& self) { return self.block_len; })
        .def_property_readonly("ecc_len", [](NTT& self) { return self.ecc_len; })
        .def_property_readonly("ecc_size", [](NTT& self) { return self.ecc_len * sizeof(uint16_t); })

        .def("rbo", &NTT::rbo, R"(Reverse Bit Order for :math:`\log_2 \text{ntt_len}` bits)", "i"_a)

        .def("ntt", cast_args(&NTT::py_ntt<simd_map_t<1>>), R"(RBO input, normal order output)", "input"_a)
        .def("ntt4", cast_args(&NTT::py_ntt<simd_map_t<4>>), R"(RBO input, normal order output)", "input"_a)
        .def("ntt8", cast_args(&NTT::py_ntt<simd_map_t<8>>), R"(RBO input, normal order output)", "input"_a)
        .def("ntt16", cast_args(&NTT::py_ntt<simd_map_t<16>>), R"(RBO input, normal order output)", "input"_a)

        .def("intt", cast_args(&NTT::py_intt<simd_map_t<1>>), R"(Normal order input, RBO output)", "input"_a)
        .def("intt4", cast_args(&NTT::py_intt<simd_map_t<4>>), R"(Normal order input, RBO output)", "input"_a)
        .def("intt8", cast_args(&NTT::py_intt<simd_map_t<8>>), R"(Normal order input, RBO output)", "input"_a)
        .def("intt16", cast_args(&NTT::py_intt<simd_map_t<16>>), R"(Normal order input, RBO output)", "input"_a)

        .def("nttr", cast_args(&NTT::py_nttr<simd_map_t<1>>), R"(Normal order input, RBO output)", "input"_a)
        .def("nttr4", cast_args(&NTT::py_nttr<simd_map_t<4>>), R"(Normal order input, RBO output)", "input"_a)
        .def("nttr8", cast_args(&NTT::py_nttr<simd_map_t<8>>), R"(Normal order input, RBO output)", "input"_a)
        .def("nttr16", cast_args(&NTT::py_nttr<simd_map_t<16>>), R"(Normal order input, RBO output)", "input"_a)

        .def("inttr", cast_args(&NTT::py_inttr<simd_map_t<1>>), R"(RBO input, Normal order output)", "input"_a)
        .def("inttr4", cast_args(&NTT::py_inttr<simd_map_t<4>>), R"(RBO input, Normal order output)", "input"_a)
        .def("inttr8", cast_args(&NTT::py_inttr<simd_map_t<8>>), R"(RBO input, Normal order output)", "input"_a)
        .def("inttr16", cast_args(&NTT::py_inttr<simd_map_t<16>>), R"(RBO input, Normal order output)", "input"_a)

        .def("poly_mul", cast_args(&NTT::py_poly_mul<simd_map_t<1>>), R"(Element wise multiplication)", "a"_a, "b"_a)
        .def("poly_mul4", cast_args(&NTT::py_poly_mul<simd_map_t<4>>), R"(Element wise multiplication x4)", "a"_a, "b"_a)
        .def("poly_mul8", cast_args(&NTT::py_poly_mul<simd_map_t<8>>), R"(Element wise multiplication x8)", "a"_a, "b"_a)
        .def("poly_mul16", cast_args(&NTT::py_poly_mul<simd_map_t<16>>), R"(Element wise multiplication x16)", "a"_a, "b"_a)

        .def("poly_div", cast_args(&NTT::py_poly_div<simd_map_t<1>>), R"(Polynomial division)", "a"_a, "b"_a)
        .def("poly_div4", cast_args(&NTT::py_poly_div<simd_map_t<4>>), R"(Polynomial division x4)", "a"_a, "b"_a)
        .def("poly_div8", cast_args(&NTT::py_poly_div<simd_map_t<8>>), R"(Polynomial division x8)", "a"_a, "b"_a)
        .def("poly_div16", cast_args(&NTT::py_poly_div<simd_map_t<16>>), R"(Polynomial division x16)", "a"_a, "b"_a)

        .def("poly_inv", cast_args(&NTT::py_poly_inv<simd_map_t<1>>), R"(Polynomial inverse)", "a"_a, "n"_a)
        .def("poly_inv4", cast_args(&NTT::py_poly_inv<simd_map_t<4>>), R"(Polynomial inverse x4)", "a"_a, "n"_a)
        .def("poly_inv8", cast_args(&NTT::py_poly_inv<simd_map_t<8>>), R"(Polynomial inverse x8)", "a"_a, "n"_a)
        .def("poly_inv16", cast_args(&NTT::py_poly_inv<simd_map_t<16>>), R"(Polynomial inverse x16)", "a"_a, "n"_a)
    ;
}
