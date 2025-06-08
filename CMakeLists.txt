cmake_minimum_required(VERSION 3.4...3.18)
project(ffrs)

add_subdirectory(pybind11)


if (MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -O3 -flto=auto)
endif()


pybind11_add_module(pyffrs "libffrs/pyffrs.cpp")
set_target_properties(pyffrs PROPERTIES OUTPUT_NAME "libffrs")
target_compile_features(pyffrs PRIVATE cxx_std_17)
target_compile_definitions(pyffrs PRIVATE VERSION_INFO="0.1.0")

target_precompile_headers(pyffrs PRIVATE
    <pybind11/pybind11.h>
    <pybind11/stl.h>
)
