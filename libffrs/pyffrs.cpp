/**************************************************************************
 * pyffrs.cpp
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

#include <optional>

#include <pybind11/pybind11.h>

#include "pygf256.hpp"
#include "pygfi16.hpp"
#include "pyrs256.hpp"
#include "pyrsi16.hpp"
#include "pycirc16.hpp"
#include "pyntt.hpp"

#include <pybind11/pybind11.h>

#include <sys/mman.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>
#include <stdexcept>


// TODO: make portable
py::memoryview py_create_buffer(size_t requested_size) {
    // round up to hugepage size
    constexpr size_t hugepage_size = 2 * 1024 * 1024;  // 2 MiB
    size_t alloc_size = ((requested_size + hugepage_size - 1) / hugepage_size) * hugepage_size;

    void *ptr = mmap(
        nullptr,
        alloc_size,
        PROT_READ | PROT_WRITE,
        MAP_PRIVATE | MAP_ANONYMOUS, // | MAP_HUGETLB,
        -1,
        0
    );

    if (ptr == MAP_FAILED) {
        throw std::runtime_error(
            std::string("mmap failed: ") + std::strerror(errno));
    }

    madvise(ptr, alloc_size, MADV_HUGEPAGE);
    mlock(ptr, alloc_size);
    memset(ptr, 0, alloc_size);

    auto *mapping = new std::pair{ptr, alloc_size};
    py::capsule capsule(
        mapping,
        "hugepage_array",
        [](PyObject *capsule) {
            auto *mapping = reinterpret_cast<std::pair<void *, size_t> *>(PyCapsule_GetPointer(capsule, "hugepage_array"));

            if (mapping) {
                munmap(mapping->first, mapping->second);
                delete mapping;
            }
        }
    );

    // return py::array_t<uint8_t>(
    //     {static_cast<py::ssize_t>(requested_size)},
    //     {static_cast<py::ssize_t>(sizeof(uint8_t))},
    //     static_cast<uint8_t *>(ptr),
    //     capsule
    // );

    Py_buffer view;
    auto result = PyBuffer_FillInfo(
        &view,
        capsule.ptr(),  // owner object
        ptr,            // data pointer
        requested_size, // size in bytes
        0,              // writable
        PyBUF_CONTIG    // flags
    );
    if (result < 0)
        throw py::error_already_set();

    return py::reinterpret_steal<py::memoryview>(PyMemoryView_FromBuffer(&view));
}


PYBIND11_MODULE(libffrs, m) {
    using namespace pybind11::literals;

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

    m.def("create_buffer", &py_create_buffer, "size"_a, R"(
        Create a memory buffer of the specified size, backed by hugepages if possible.

        The buffer is returned as a memoryview object that can be used in Python.
    )");

    m.def("set_logger", [](py::object&& logger) { PyLogger::set_logger(logger); },
        "logger"_a, R"(
        Logger object to be used by C++ library or ``None`` to disable logging.
    )");

    m.doc() = R"(
        FFRS - Fairly Fast & Flexible Reed-Solomon coding
    )";

    PyGF256::register_class(m);
    PyRS256::register_class(m);
    PyGFi16::register_class(m);
    PyRSi16::register_class(m);
    PyCIRC16::register_class(m);
    NTT::register_class(m);
}
