/**************************************************************************
 * util.hpp
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

#pragma once

#include <pybind11/pybind11.h>


namespace detail {

/**
 * Adaptor to automatically request buffer_info and cast pointer to data
 */
template<typename T, bool writeable>
class buffer_adpator {
public:
    pybind11::buffer buffer;
    pybind11::buffer_info info;
    size_t size = 0;
    T *data = nullptr;

    inline buffer_adpator(const pybind11::buffer& o):
        buffer(std::move(o)),
        info(buffer.request(false)),
        size(size_t(info.size) / sizeof(T)),
        data(reinterpret_cast<T *>(info.ptr))
    {
        if (info.size % sizeof(T) != 0)
            throw pybind11::value_error("Invalid buffer size: " + std::to_string(info.size));
    }

    inline T& operator[](size_t idx) {
        return data[idx];
    }
};

}

template<typename T>
using buffer_ro = detail::buffer_adpator<const T, false>;
template<typename T>
using buffer_rw = detail::buffer_adpator<T, true>;


namespace detail {

template<typename T, typename = std::remove_reference_t<T>>
struct cast_arg { using type = T; };

template<typename T, typename U>
struct cast_arg<T, buffer_ro<U>> { using type = pybind11::buffer&&; };

template<typename T, typename U>
struct cast_arg<T, buffer_rw<U>> { using type = pybind11::buffer&&; };

template<typename T>
using cast_arg_t = typename cast_arg<T>::type;

}

template <typename Return, typename Class, typename...Args>
constexpr auto cast_args(Return (Class::*method)(Args...)) {
    return [method](Class& self, typename detail::cast_arg<Args>::type...args) {
        return (self.*method)(std::move(args)...);
    };
}

#define _PY_ASSERT_MSG(c, text) do { if (!(c)) throw std::runtime_error(__FILE__ ":" + std::to_string(__LINE__) + ": Assertion failed: " #c ", " + text); } while(0)
#define _PY_ASSERT(c)           do { if (!(c)) throw std::runtime_error(__FILE__ ":" + std::to_string(__LINE__) + ": Assertion failed: " #c); } while(0)
#define _PY_ASSERT_ARG(arg1, arg2, arg3, ...) arg3
#define _PY_ASSERT_GET(...) _PY_ASSERT_ARG(__VA_ARGS__, _PY_ASSERT_MSG, _PY_ASSERT)
#define py_assert(...) _PY_ASSERT_GET(__VA_ARGS__)(__VA_ARGS__)