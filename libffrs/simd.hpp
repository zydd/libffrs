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


namespace vec {

template<typename GFT>
simd_mask_t is_zero(GFT const& vec);

template<typename GFT>
void copy_n_masked(const GFT *const src, size_t n, GFT *const dst, GFT const& condition);

// template<typename GFT>
// void assign_masked(GFT& vec, GFT const& value, GFT const& condition);


template<typename GFT>
inline bool any(GFT const& a) {
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


template<typename GFT>
inline ::GFT min(GFT const& a) {
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


template<typename GFT>
inline GFT min(GFT const& a, GFT const& b) {
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


template<typename GFT>
inline ::GFT max(GFT const& a) {
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


template<typename GFT>
inline GFT max(GFT const& a, GFT const& b) {
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

template<typename GFT>
inline GFT gather(const GFT *const src, GFT const& i) {
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


template<typename GFT>
void copy_n(const GFT *const src, GFT const& src_offset, GFT n, GFT *const dst, GFT const& dst_offset) {
    if constexpr (std::is_integral_v<GFT>) {
        std::copy_n(&src[src_offset], n, &dst[dst_offset]);
    } else {
        auto max_n = vec::max(n);
        for (size_t i = 0; i < max_n; i++) {
            for (size_t j = 0; j < sizeof(GFT) / sizeof(::GFT); j++)
                if (i < n[j])
                    dst[i + dst_offset[j]][j] = src[i + src_offset[j]][j];
        }
    }
}


template<typename GFT>
void fill_n(GFT *const dst, GFT const& dst_offset, GFT n, GFT const& value) {
    if constexpr (std::is_integral_v<GFT>) {
        std::fill_n(&dst[dst_offset], n, value);
    } else {
        auto max_n = vec::max(n);
        for (size_t i = 0; i < max_n; i++) {
            for (size_t j = 0; j < sizeof(GFT) / sizeof(::GFT); j++)
                if (i < n[j])
                    dst[i + dst_offset[j]][j] = value[j];
        }
    }
}

template<typename Src, typename Dst>
inline void copy_transposed(const Src src[], size_t src_cols, Dst dst[], size_t dst_cols) {
    for (size_t j = 0; j < src_cols; ++j)
        for (size_t i = 0; i < dst_cols; ++i)
            dst[i + j * dst_cols] = src[i * src_cols + j];
}


template<typename Src, typename Dst, bool = std::is_integral_v<Src>, bool = std::is_integral_v<Dst>>
inline void copy_transposed(const Src src[], size_t src_stride, size_t src_cols, Dst dst[], size_t dst_stride, size_t dst_cols) {
    for (size_t j = 0; j < src_cols; ++j)
        for (size_t i = 0; i < dst_cols; ++i)
            dst[i + j * dst_stride] = src[i * src_stride + j];
}


template<bool Prefetch = true, size_t PrefetchDistance = 0, typename Src, typename Dst>
inline void copy_stride(const Src *src, size_t src_stride, Dst *dst, size_t dst_stride, size_t cols, size_t count) {
    size_t src_i = 0;
    size_t dst_i = 0;
    for (size_t i = 0; i < count; ++i) {
        if constexpr (Prefetch) {
            auto next = src_i + src_stride;
            __builtin_prefetch(&src[next]);

            if constexpr (PrefetchDistance)
                __builtin_prefetch(&src[next + PrefetchDistance]);
        }
        std::copy_n(&src[src_i], cols, &dst[dst_i]);
        // std::memcpy(&dst[dst_i], &src[src_i], cols * sizeof(T));
        src_i += src_stride;
        dst_i += dst_stride;
    }

    // for (size_t i = 0; i < count; ++i) {
    //     // std::copy_n(src + i * src_stride, n, dst + i * dst_stride);
    //     std::memcpy(dst + i * dst_stride, src + i * src_stride, cols * sizeof(T));
    // }
}

}


namespace vec::poly {
/**
 * Poly degree + 1
 */
template<typename GFT>
inline GFT len(const GFT *const r, size_t r_len_max) {
    // check_le_ecc(r_len_max);

    GFT len = GFT{} + ::GFT(r_len_max);

    if constexpr (std::is_integral_v<GFT>) {
        while (r_len_max > 0 && vec::is_zero(r[r_len_max - 1]))
            --r_len_max;
        return r_len_max;
    } else {
        GFT searching = (len > 0);
        while (r_len_max > 0) {
            searching &= (r[r_len_max - 1] == 0);
            len += searching;
            searching &= (len > 0);
            --r_len_max;
        }
    }

    // check_le_ecc(len);
    return len;
}

}
