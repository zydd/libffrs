/**************************************************************************
 * detail.hpp
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

#pragma once

#include <cstdint>

// #ifdef _MSC_VER
// #define _noinline __declspec(noinline)
// #else
// #define _noinline __attribute__((noinline))
// #endif

namespace ffrs {

using size_t = std::size_t;

namespace detail {

template<typename T, unsigned Bits>
struct bit_array {
    static_assert(Bits <= sizeof(T) * 8);

    uint8_t *ptr;

    struct bit_field {
        uint8_t *ptr;
        uint8_t offset;

        void operator=(T const& in) {
            auto ptr = this->ptr;
            auto offset = this->offset;

            for (unsigned i = 0; i < Bits;) {
                uint8_t val;
                uint8_t mask = 0xff >> offset;

                uint8_t rem_bits = Bits - i;
                uint8_t avail_bits = 8 - offset;

                if (rem_bits < avail_bits) {
                    val = in << (avail_bits - rem_bits);
                    mask &= 0xff << (avail_bits - rem_bits);
                } else {
                    val = in >> (rem_bits - avail_bits);
                }

                *ptr = (*ptr & ~mask) | (val & mask);
                ++ptr;
                i += avail_bits;
                offset = 0;
            }
        }

        operator T() const {
            auto ptr = this->ptr;
            auto offset = this->offset;
            T ret = 0;

            for (unsigned i = 0; i < Bits;) {
                uint8_t mask = 0xff >> offset;

                uint8_t rem_bits = Bits - i;
                uint8_t avail_bits = 8 - offset;

                if (rem_bits < avail_bits) {
                    mask &= 0xff << (avail_bits - rem_bits);

                    ret <<= rem_bits;
                    ret |= (*ptr++ & mask) >> (avail_bits - rem_bits);
                } else {
                    ret <<= avail_bits;
                    ret |= (*ptr++ & mask);
                }

                i += avail_bits;
                offset = 0;
            }
            return ret;
        }
    };

    bit_field operator[](int i) {
        auto bit_i = i * Bits;
        return bit_field{ptr + (bit_i >> 3), uint8_t(bit_i & 0x07)};
    }
};


// Curiously Recurring Template Pattern
template<typename...Fs>
class CRTP : public Fs... {
public:
    template<typename F>
    static constexpr CRTP const& cast(F _this) {
        return static_cast<CRTP const&>(*_this);
    }

    template<typename F>
    static constexpr CRTP& cast_mut(F _this) {
        return static_cast<CRTP&>(*_this);
    }

    template<typename...F0>
    inline CRTP(F0&&... f0):
        F0(std::move(f0))...
    { }
};

template<typename T>
constexpr bool is_power_of_two(T a) {
    return a && (a & (a - 1)) == 0;
}


template<typename T>
constexpr T ilog2_floor(T a) {
    T r = 0;
    while (a >>= 1)
        r += 1;
    return r;
}


constexpr size_t ipow(size_t a, size_t b) {
    size_t r = 1;
    for (size_t i = 0; i < b; ++i)
        r *= a;
    return r;
}

constexpr size_t align_size(size_t size, size_t alignment) {
    return size + (alignment - (size % alignment)) % alignment;
}

}
}
