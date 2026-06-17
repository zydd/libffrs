/**************************************************************************
 * simd.hpp
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

# pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

#include "pygfi16.hpp"
#include "rsi16v.hpp"
#include "ntt.hpp"
#include "pylogging.hpp"

#include <pybind11/stl.h>


typedef GFT GFTx4 __attribute__((vector_size(4 * sizeof(GFT))));
typedef GFT GFTx8 __attribute__((vector_size(8 * sizeof(GFT))));
typedef GFT GFTx16 __attribute__((vector_size(16 * sizeof(GFT))));


template<int N>
using simd_map_t =
    std::conditional_t<N == 1, GFT,
    std::conditional_t<N == 4, GFTx4,
    std::conditional_t<N == 8, GFTx8,
    std::conditional_t<N == 16, GFTx16,
    void>>>>;

using simd_mask_t = int;


template<typename GFT>
struct vec {
    static void assign_masked(GFT& vec, GFT const& value, GFT const& condition);
    static void copy_n_masked(const GFT *const src, size_t n, GFT *const dst, GFT const& condition);

    static void copy_n(const GFT *const src, GFT const& src_offset, GFT n, GFT *const dst, GFT const& dst_offset) {
        if constexpr (std::is_integral_v<GFT>) {
            std::copy_n(&src[src_offset], n, &dst[dst_offset]);
        } else {
            auto max_n = max(n);
            for (size_t i = 0; i < max_n; i++) {
                for (size_t j = 0; j < sizeof(GFT) / sizeof(::GFT); j++)
                    if (i < n[j])
                        dst[i + dst_offset[j]][j] = src[i + src_offset[j]][j];
            }
        }
    }

    static void fill_n(GFT *const dst, GFT const& dst_offset, GFT n, GFT const& value) {
        if constexpr (std::is_integral_v<GFT>) {
            std::fill_n(&dst[dst_offset], n, value);
        } else {
            auto max_n = max(n);
            for (size_t i = 0; i < max_n; i++) {
                for (size_t j = 0; j < sizeof(GFT) / sizeof(::GFT); j++)
                    if (i < n[j])
                        dst[i + dst_offset[j]][j] = value[j];
            }
        }
    }

    static inline bool any(GFT const& a) {
        if constexpr (std::is_integral_v<GFT>) {
            return a != 0;
        } else {
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                if (a[i] != 0)
                    return true;
            }
            return false;
        }
    }

    static inline ::GFT min(GFT const& a) {
        if constexpr (std::is_integral_v<GFT>) {
            return a;
        } else {
            auto min = ~::GFT{0};
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                if (a[i] < min)
                    min = a[i];
            }
            return min;
        }
    }

    static inline ::GFT max(GFT const& a) {
        // *std::max_element(reinterpret_cast<::GFT *>(&root_count), reinterpret_cast<::GFT *>(&root_count + 1));

        if constexpr (std::is_integral_v<GFT>) {
            return a;
        } else {
            ::GFT max = 0;
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                if (a[i] > max)
                    max = a[i];
            }
            return max;
        }
    }

    static inline GFT min(GFT const& a, GFT const& b) {
        if constexpr (std::is_integral_v<GFT>) {
            return std::min(a, b);
        } else {
            GFT min = GFT{0};
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                min[i] = std::min(a[i], b[i]);
            }
            return min;
        }
    }

    static inline GFT max(GFT const& a, GFT const& b) {
        if constexpr (std::is_integral_v<GFT>) {
            return std::max(a, b);
        } else {
            GFT max = GFT{0};
            for (size_t i = 0; i < sizeof(GFT) / sizeof(::GFT); ++i) {
                max[i] = std::max(a[i], b[i]);
            }
            return max;
        }
    }

    static inline GFT gather(const GFT *const src, GFT const& i) {
        if constexpr (std::is_integral_v<GFT>) {
            return src[i];
        } else {
            GFT res = GFT{0};
            for (size_t j = 0; j < sizeof(GFT) / sizeof(::GFT); ++j) {
                res[j] = src[i[j]][j];
            }
            return res;
        }
    }

    static simd_mask_t is_zero(GFT const& vec);
};
