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

#include "pygf256.hpp"
#include "pygfi32.hpp"
#include "pyrs256.hpp"
#include "pyrsi16.hpp"


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

    PyGF256::register_class(m);
    PyRS256::register_class(m);
    PyGFi32::register_class(m);
    PyRSi16::register_class(m);
}
