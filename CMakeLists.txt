cmake_minimum_required(VERSION 3.4...3.18)
project(ffrs)

add_subdirectory(pybind11)

pybind11_add_module(pyffrs "src/pyffrs.cpp")
set_target_properties(pyffrs PROPERTIES OUTPUT_NAME "ffrs")
target_include_directories(pyffrs PRIVATE "./include")
target_compile_features(pyffrs PRIVATE cxx_std_17)
target_compile_definitions(pyffrs PRIVATE VERSION_INFO="0.1.0")
