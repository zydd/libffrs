/**************************************************************************
 * detail.hpp
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

namespace ffrs {
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

template<typename T, typename...Fs>
class CRTP : public Fs... {
public:
    template<typename F>
    static constexpr T const& cast(F _this) {
        return static_cast<T const&>(*_this);
    }

protected:
    // base class constructors are always called first
    // define an `init` method to make sure it's called after GF is initialized
    // this method needs to be called from the subclasses constructors
    inline void init() { init_mixins<Fs...>(); }

private:
    template <typename Fs0>
    auto try_init_mixin(int) -> decltype(Fs0::init()) { return this->Fs0::init(); }

    template <typename Fs0>
    auto try_init_mixin(long) -> void { }

    template<typename Fs0, typename...Fss>
    inline void init_mixins() {
        try_init_mixin<Fs0>(0);

        if constexpr (sizeof...(Fss) > 0)
            init_mixins<Fss...>();
    }
};


inline constexpr unsigned ilog2_floor(unsigned a) {
    unsigned r = 0;
    while (a >>= 1)
        r += 1;
    return r;
}


inline constexpr unsigned ipow(unsigned a, unsigned b) {
    unsigned r = 1;
    for (unsigned i = 0; i < b; ++i)
        r *= a;
    return r;
}

}
}
