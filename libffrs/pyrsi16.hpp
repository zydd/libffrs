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

    inline rs_encode_ntt() {
        auto& gf = RS::cast(this).gf;
        for (unsigned i = 0; i < MaxFieldBits; ++i) {
            GFT n = 1 << i;
            if (n >= gf.field_elements)
                break;

            GFT root = gf.exp(gf.div(gf.log(1), n));
            if (gf.pow(root, n) != 1)
                break;

            _nth_roots_of_unity[i] = root;
        }
    }

    inline GFT get_root_of_unity(size_t n) {
        return  _nth_roots_of_unity[__builtin_ctzl(n)];
    }

    // inline void encode(const uint8_t input[], size_t input_size, uint8_t output[]) {
    // }

protected:
    GFT _nth_roots_of_unity[MaxFieldBits] = {0};
};

template<typename GF>
using RSi16 = ffrs::RS<GF, ffrs::rs_data,
    rs_encode_ntt
    >;


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
        RSi16<PyGFi32>(rs_data(PyGFi32(65537, primitive), ecc_len))
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
        auto root = get_root_of_unity(default_block_size / sizeof(uint16_t));
        size_t msg_size = buf.size - ecc_len / sizeof(uint16_t);

        std::vector<uint16_t> temp(buf.size);
        gf.ntt(root, &buf[0], buf.size, &temp[0], buf.size);
        std::copy_n(&temp[0], ecc_len / sizeof(uint16_t), &buf[msg_size]);
    }

    inline py::bytearray py_encode_blocks(buffer_ro<uint16_t> buf) {
        size_t output_block_size = default_block_size / sizeof(uint16_t);
        size_t ecc_len = this->ecc_len / sizeof(uint16_t);
        size_t input_block_size = output_block_size - ecc_len;
        auto root = get_root_of_unity(output_block_size);

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
            gf.ntt(root, &buf.data[block * input_block_size], input_block_size, &output_block[input_block_size], output_block_size);
            std::copy_n(&buf[block * input_block_size], input_block_size, &output_block[0]);
        }

        if (input_remainder > 0) {
            auto output_block = &output_data[output_size - input_remainder - ecc_len];
            gf.ntt(root, &buf[buf.size - input_remainder], input_remainder, &output_block[input_remainder], output_block_size);
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

