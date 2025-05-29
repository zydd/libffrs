/**************************************************************************
 * pyffrs.cpp
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

#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ntt.hpp"
#include "reed_solomon.hpp"

#include "util.hpp"

#include "pygf256.hpp"
#include "pyrs256.hpp"
#include "pygfi32.hpp"

namespace py = pybind11;
using namespace pybind11::literals;



template<typename GF>
using NTT256x8 = ffrs::NTT<GF,
    ffrs::ntt_data,
    // ffrs::ntt_eval
    ffrs::ntt_eval_lut<256>::type
    >;



class PyNTT256x8 : public NTT256x8<PyGF256> {
public:
    inline PyNTT256x8(uint8_t primitive, uint16_t poly1, std::optional<uint8_t> root):
        NTT256x8<PyGF256>(ntt_data(PyGF256(primitive, poly1), root.value_or(primitive)))
    { }

    inline py::bytearray py_ntt8(buffer_ro<uint8_t> buf) {
        uint64_t res = ntt(buf.data, buf.size);
        return py::bytearray(reinterpret_cast<const char *>(&res), std::min(buf.size, sizeof(uint64_t)));
    }
};



PYBIND11_MODULE(libffrs, m) {
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

    py::class_<PyNTT256x8>(m, "NTT256x8")
        .def_property_readonly("gf", [](PyNTT256x8& self) -> auto const& { return self.gf; })

        .def(py::init<uint8_t, uint16_t, std::optional<uint8_t>>(), R"(
            Instantiate type NTT computation over :math:`GF(p^n)/P`

            Args:
                primitive : :math:`a` -- primitive value used to generate the field
                polynomial : :math:`P` -- irreducible polynomial used to generate the field
                root : :math:`\alpha` -- root of unity for NTT computation
            )",
            "primitive"_a = 2, "poly1"_a = 0x11d, "root"_a = py::none()
        )
        .def("ntt8", cast_args(&PyNTT256x8::py_ntt8), R"(
            Number Theoretical Transform
        )", "buffer"_a)
        .doc() = R"(
            NTT for :math:`GF(2^8)`
        )";

    PyGF256::register_class(m);
    PyRS256::register_class(m);
    PyGFi32::register_class(m);
}
