/**************************************************************************
 * pyrsi16.hpp
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

# pragma once

#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "reed_solomon.hpp"
#include "util.hpp"
#include "pygfi32.hpp"

namespace py = pybind11;


template<typename GF, typename RS, size_t MaxFieldBits = 16>
class rs_encode_ntt {
public:
    using GFT = typename GF::GFT;

    inline rs_encode_ntt(size_t block_size) {
        auto& rs = RS::cast(this);
        py_assert(rs.gf.prime > 0);

        for (unsigned i = 0; i < MaxFieldBits; ++i) {
            GFT n = 1 << i;
            if (n >= rs.gf.field_elements)
                break;

            GFT root = rs.gf.exp(rs.gf.div(rs.gf.log(1), n));
            if (rs.gf.pow(root, n) != 1)
                break;

            _nth_roots_of_unity[i] = root;
        }

        GFT r = get_root_of_unity(block_size / sizeof(uint16_t));
        for (size_t i = 0; i < MaxFieldBits; ++i) {
            _roots[i] = r;
            r = rs.gf.pow(r, 2);
        }
    }

    inline GFT get_root_of_unity(size_t n) const {
        return  _nth_roots_of_unity[__builtin_ctzl(n)];
    }

    template<typename T, typename U>
    inline void encode(const T input[], size_t input_size, U output[], size_t output_size) {
        py_assert(output_size >= input_size);

        std::copy_n(input, input_size, output);
        // auto& gf = RS::cast(this).gf;
        // gf.copy_rbo(input, input_size, output, output_size);
        ct_butterfly(output, output_size);
    }

protected:
    GFT _nth_roots_of_unity[MaxFieldBits] = {0};
    GFT _roots[65536] = {0};

    template<typename T>
    inline void ct_butterfly(T data[], size_t output_size) const {
        GFT root = get_root_of_unity(output_size);
        auto& gf = RS::cast(this).gf;
        for (size_t stride = 1, exp_f = output_size >> 1; stride < output_size; stride <<= 1, exp_f >>= 1) {
            for (size_t start = 0; start < output_size; start += stride * 2) {
                // For each pair of the CT butterfly operation.
                for (size_t i = start; i < start + stride; ++i) {
                    // j = i - start
                    GFT w = gf.pow(root, exp_f * (i - start));
                    // GFT w = _roots[exp_f * (i - start)];

                    // Cooley-Tukey butterfly
                    T a = data[i];
                    T b = data[i+stride];
                    GFT m = gf.mul(w, b);
                    data[i] = gf.add(a, m);
                    data[i+stride] = gf.sub(a, m);
                }
            }
        }
    }

    template<typename T, typename U>
    inline void copy_rbo(const T input[], size_t input_size, U output[], size_t output_size) const {
        output[0] = input[0];

        if (output_size <= 2)
            return;

        size_t l, rev = 0;
        for (size_t i = 1; i < input_size; ++i) {
            for (l = output_size >> 1; rev + l >= output_size; l >>= 1);
            rev = (rev & (l - 1)) + l;
            // rev = rbo16(i) >> (16 - nbits);
            output[rev] = input[i];
        }
    }

    uint16_t rbo16(uint16_t b) const {
        b = ((b >> 1) & 0x5555) | ((b & 0x5555) << 1);
        b = ((b >> 2) & 0x3333) | ((b & 0x3333) << 2);
        b = ((b >> 4) & 0x0f0f) | ((b & 0x0f0f) << 4);
        b = ((b >> 8)         ) | ((b         ) << 8);
        return b;
    }
};

template<typename GF>
using RSi16 = ffrs::RS<GF, rs_encode_ntt, ffrs::rs_data>;


class PyRSi16 : public RSi16<PyGFi32> {
public:
    size_t default_block_size = 255;

    inline PyRSi16(
            std::optional<uint16_t> block_size,
            std::optional<uint16_t> message_len,
            std::optional<uint16_t> ecc_len,
            uint16_t primitive):
        PyRSi16(_get_ecc_len(block_size, message_len, ecc_len),
                block_size.value_or(255),
                primitive)
    {
        if (ecc_len && message_len) {
            if (block_size && *message_len + *ecc_len != *block_size) {
                throw py::value_error("block_len must be equal to message_len + ecc_len");
            } else if (!block_size) {
                set_default_block_size(size_t(*message_len) + size_t(*ecc_len));
            }
        }
    }

    inline PyRSi16(uint16_t ecc_len, size_t block_size, uint16_t primitive):
        RSi16<PyGFi32>(rs_encode_ntt(block_size), rs_data(PyGFi32(65537, primitive), ecc_len))
    {
        set_default_block_size(block_size);
    }

    inline void set_default_block_size(size_t block_size) {
        if (block_size <= ecc_len)
            throw py::value_error("block_len must be greater than ecc_len");

        if (block_size > 65536)
            throw py::value_error("block_len must be <= 65536");

        default_block_size = block_size;
    }

    inline void py_encode(buffer_rw<uint16_t> buf) {
        size_t msg_size = buf.size - ecc_len / sizeof(uint16_t);

        std::vector<uint16_t> temp(buf.size);
        encode(&buf[0], buf.size, &temp[0], buf.size);
        std::copy_n(&temp[0], ecc_len / sizeof(uint16_t), &buf[msg_size]);
    }

    inline py::bytearray py_encode_blocks(buffer_ro<uint16_t> buf) {
        size_t output_block_size = default_block_size / sizeof(uint16_t);
        size_t ecc_len = this->ecc_len / sizeof(uint16_t);
        size_t input_block_size = output_block_size - ecc_len;

        if (buf.size == 0 || output_block_size == 0 || output_block_size <= ecc_len)
            return {};

        size_t full_blocks = buf.size / input_block_size;

        size_t output_size = full_blocks * (input_block_size + ecc_len);

        // Last block will be smaller if input size is not divisible by input_block_size
        size_t input_remainder = buf.size - full_blocks * input_block_size;
        if (input_remainder > 0)
            output_size += input_remainder + ecc_len;

        // Allocate extra space to compute NTT
        auto output = py::bytearray(nullptr, (output_size + input_block_size) * sizeof(uint16_t));
        auto output_data = reinterpret_cast<uint16_t *>(PyByteArray_AsString(output.ptr()));
        std::fill_n(output_data, output_size + input_block_size, 0);

        for (size_t block = 0; block < buf.size / input_block_size; ++block) {
            auto output_block = &output_data[block * output_block_size];
            encode(&buf.data[block * input_block_size], input_block_size, &output_block[input_block_size], output_block_size);
            std::copy_n(&buf[block * input_block_size], input_block_size, &output_block[0]);
        }

        if (input_remainder > 0) {
            auto output_block = &output_data[output_size - input_remainder - ecc_len];
            encode(&buf[buf.size - input_remainder], input_remainder, &output_block[input_remainder], output_block_size);
            std::copy_n(&buf[buf.size - input_remainder], input_remainder, output_block);
        }

        // Remove extra space
        PyByteArray_Resize(output.ptr(), output_size * sizeof(uint16_t));
        return output;
    }

    // inline bool py_decode(buffer_rw<uint8_t> buf) {
    //     size_t msg_size = buf.size - ecc_len;
    //     return decode(buf, msg_size, &buf[msg_size]);
    // }

    // inline py::bytearray py_synds(buffer_ro<uint8_t> buf) {
    //     if (buf.size < ecc_len)
    //         return {};

    //     size_t msg_size = buf.size - ecc_len;
    //     synds_array_t synds_arr;

    //     synds(buf, msg_size, &buf[msg_size], synds_arr);

    //     return py::bytearray(reinterpret_cast<const char *>(synds_arr), ecc_len);
    // }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyRSi16>(m, "RSi16")
            .def_property_readonly("ecc_len", [](PyRSi16& self) { return self.ecc_len; })

            .def_property("block_len",
                [](PyRSi16& self) { return self.default_block_size; },
                &PyRSi16::set_default_block_size)

            .def_property_readonly("message_len",
                [](PyRSi16& self) { return self.default_block_size - self.ecc_len; })

            .def_property_readonly("gf", [](PyRSi16& self) -> auto const& { return self.gf; })

            .def_property_readonly("roots_of_unity", [](PyRSi16& self) {
                return std::vector<uint16_t>(std::begin(self._nth_roots_of_unity), std::end(self._nth_roots_of_unity)); })

            .def(py::init<std::optional<uint16_t>, std::optional<uint16_t>, std::optional<uint16_t>, uint16_t>(), R"(
                Instantiate a Reed-Solomon encoder with the given configuration)",
                "block_len"_a = py::none(), "message_len"_a = py::none(), "ecc_len"_a = py::none(), "primitive"_a = 3)

            .def("__sizeof__", [](PyRSi16& self) { return sizeof(self); })

            .def("encode", cast_args(&PyRSi16::py_encode),
                R"(Systematic encode)",
                "buffer"_a)

            .def("encode_blocks", cast_args(&PyRSi16::py_encode_blocks), R"(
                Encode the input buffer in blocks, storing the parity bytes right next to their corresponding block
                )",
                "buffer"_a)

            // .def("decode", cast_args(&PyRSi16::py_decode),
            //     R"(Systematic decode)",
            //     "buffer"_a)

            // .def("_synds", cast_args(&PyRSi16::py_synds),
            //     R"(Compute syndromes)",
            //     "buffer"_a)

            .doc() = R"(Reed-Solomon coding over :math:`GF(2^8)`)";
    }

private:
    inline uint16_t _get_ecc_len(
            std::optional<uint16_t> block_len,
            std::optional<uint16_t> message_len,
            std::optional<uint16_t> ecc_len) {

        uint16_t res;
        if (ecc_len) {
            res = *ecc_len;
        } else if (message_len && block_len) {
            if (*message_len >= *block_len) {
                throw py::value_error("block_len must be greater than message_len");
            }
            res = block_len.value_or(default_block_size) - *message_len;
        } else {
            throw py::value_error("Must specify either (block_len, message_len) or ecc_len");
        }

        if (res % 2 != 0)
            throw py::value_error("ecc_len must be a multiple of 2: " + std::to_string(res));

        return res;
    }
};

