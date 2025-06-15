cmake_minimum_required(VERSION 3.4...3.18)
project(ffrs)

add_subdirectory(pybind11)


if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O3 -flto=auto)
endif()
add_compile_options(-Wall -Wextra -Wpedantic)


pybind11_add_module(pyffrs
    "libffrs/pyffrs.cpp"
    "libffrs/rsi16md_sse2.cpp"
    "libffrs/rsi16md_avx2.cpp"
    "libffrs/rsi16md_avx512.cpp"
)

set_target_properties(pyffrs PROPERTIES OUTPUT_NAME "libffrs")
target_compile_features(pyffrs PRIVATE cxx_std_17)
target_compile_definitions(pyffrs PRIVATE VERSION_INFO="0.1.0")

set_source_files_properties("libffrs/rsi16md_sse2.cpp" PROPERTIES COMPILE_FLAGS "-msse -msse2")
set_source_files_properties("libffrs/rsi16md_avx2.cpp" PROPERTIES COMPILE_FLAGS "-mavx2")
set_source_files_properties("libffrs/rsi16md_avx512.cpp" PROPERTIES COMPILE_FLAGS "-mavx512f")
