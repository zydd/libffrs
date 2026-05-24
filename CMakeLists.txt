cmake_minimum_required(VERSION 3.4...3.18)
project(ffrs)

# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

add_subdirectory(pybind11)

pybind11_add_module(pyffrs
    "libffrs/pyffrs.cpp"
    "libffrs/rsi16v.cpp"
    "libffrs/rsi16v_sse2.cpp"
    "libffrs/rsi16v_avx2.cpp"
    "libffrs/rsi16v_avx512.cpp"
)

set_source_files_properties("libffrs/rsi16v_sse2.cpp" PROPERTIES COMPILE_FLAGS "-msse -msse2")
set_source_files_properties("libffrs/rsi16v_avx2.cpp" PROPERTIES COMPILE_FLAGS "-mavx2")
# TODO: create variant with only -mavx512f
set_source_files_properties("libffrs/rsi16v_avx512.cpp" PROPERTIES COMPILE_FLAGS "-mavx512f -mavx512dq")


set_target_properties(pyffrs PROPERTIES OUTPUT_NAME "libffrs")
target_compile_features(pyffrs PRIVATE cxx_std_20)
target_compile_definitions(pyffrs PRIVATE VERSION_INFO="0.1.0")
target_compile_options(pyffrs PRIVATE -Wall -Wextra -Wpedantic)


add_custom_command(TARGET pyffrs POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E env
        PYTHONPATH=$<TARGET_FILE_DIR:pyffrs>
        python3 ${CMAKE_SOURCE_DIR}/scripts/gen-stubs.py "${CMAKE_SOURCE_DIR}/libffrs/"
    VERBATIM
)


if(SANITIZE_ADDRESS STREQUAL "ON")
    # To run with -fsanitize=address
    # LD_PRELOAD=$(clang -print-file-name=libclang_rt.asan-x86_64.so)
    # LD_PRELOAD=$(gcc -print-file-name=libasan.so)
    # ASAN_OPTIONS=detect_leaks=0

    target_link_options(pyffrs PRIVATE
        -fno-omit-frame-pointer
        -fsanitize=address
    )
    target_link_options(pyffrs PRIVATE
        -fno-omit-frame-pointer
        -fsanitize=address
    )
endif()


if(SANITIZE_UNDEFINED STREQUAL "ON")
    # To run with -fsanitize=undefined
    # LD_PRELOAD=$(clang -print-file-name=libclang_rt.ubsan-x86_64.so)
    # LD_PRELOAD=$(gcc -print-file-name=libubsan.so)

    target_link_options(pyffrs PRIVATE
        -fsanitize=undefined
    )
    target_link_options(pyffrs PRIVATE
        -fsanitize=undefined
    )
endif()


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(pyffrs PRIVATE
        -g
        -O0

        # -fhardened  # GCC 14
        # https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html
        # equivalent to:
        -D_FORTIFY_SOURCE=3
        -D_GLIBCXX_ASSERTIONS
        -ftrivial-auto-var-init=zero
        -fstack-protector-strong
        -fstack-clash-protection
        -fcf-protection=full

        -ftrivial-auto-var-init=pattern
        -D_GLIBCXX_DEBUG
        -D_GLIBCXX_DEBUG_PEDANTIC
        # -save-temps -fverbose-asm
    )
    target_link_options(pyffrs PRIVATE
        -Wl,-z,relro,-z,now
        -fstack-protector-strong
        -fstack-clash-protection
        -fcf-protection=full
    )
else()
    target_compile_options(pyffrs PRIVATE
        -O3
        -flto=auto
        # -funroll-loops
        # -save-temps -fverbose-asm
    )
endif()
