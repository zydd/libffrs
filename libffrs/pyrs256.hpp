/**************************************************************************
 * pyrs256.hpp
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

#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "reed_solomon.hpp"
#include "util.hpp"
#include "pygf256.hpp"

namespace py = pybind11;


template<typename GF>
using RS256 = ffrs::RS<GF, ffrs::rs_data,
    ffrs::rs_generator<255>::type,

    // ffrs::rs_encode_basic,
    // ffrs::rs_encode_basic_v2,
    // ffrs::rs_encode_lut_pw2<256>::type,
    ffrs::rs_encode_slice_pw2_dispatch<255>::type,

    ffrs::rs_ntt_data,
    // // ffrs::ntt_eval,
    // ffrs::ntt_eval_lut<256>::type,
    // ffrs::rs_encode_ntt,

    // ffrs::rs_synds_basic<256>::type,
    ffrs::rs_synds_lut_pw2<uint32_t, 255>::type,

    // ffrs::rs_roots_eval_basic,
    // ffrs::rs_roots_eval_uint8_chien,
    ffrs::rs_roots_eval_lut_pw2<uint64_t>::type,

    ffrs::rs_decode
    >;


class PyRS256 : public RS256<PyGF256> {
public:
    size_t default_block_size = 255;

    inline PyRS256(
            std::optional<uint8_t> block_size,
            std::optional<uint8_t> message_len,
            std::optional<uint8_t> ecc_len,
            uint8_t primitive,
            uint16_t polynomial):
        PyRS256(_get_ecc_len(block_size, message_len, ecc_len),
                block_size.value_or(255),
                primitive, polynomial)
    {
        if (ecc_len && message_len) {
            if (block_size && *message_len + *ecc_len != *block_size) {
                throw py::value_error("block_len must be equal to message_len + ecc_len");
            } else if (!block_size) {
                set_default_block_size(size_t(*message_len) + size_t(*ecc_len));
            }
        }
    }

    inline PyRS256(uint8_t ecc_len, size_t block_size, uint8_t primitive, uint16_t polynomial):
        RS256<PyGF256>(rs_data(PyGF256(primitive, polynomial), ecc_len), rs_ntt_data(primitive))
    {
        set_default_block_size(block_size);
    }

    inline void set_default_block_size(size_t block_size) {
        if (block_size <= ecc_len)
            throw py::value_error("block_len must be greater than ecc_len");

        if (block_size > 255)
            throw py::value_error("block_len must be <= 255");

        default_block_size = block_size;
    }

    inline void py_encode(buffer_rw<uint8_t> buf) {
        size_t msg_size = buf.size - ecc_len;
        encode(buf.data, msg_size, &buf.data[msg_size]);
    }

    inline py::bytearray py_encode_blocks(buffer_ro<uint8_t> buf, std::optional<size_t> block_size) {
        size_t output_block_size = block_size.value_or(default_block_size);
        size_t input_block_size = output_block_size - ecc_len;

        if (buf.size == 0 || output_block_size == 0 || output_block_size <= ecc_len)
            return {};

        size_t full_blocks = buf.size / input_block_size;

        size_t output_size = full_blocks * (input_block_size + this->ecc_len);

        // Last block will be smaller if input size is not divisible by input_block_size
        size_t input_remainder = buf.size - full_blocks * input_block_size;
        if (input_remainder > 0)
            output_size += input_remainder + this->ecc_len;

        auto output = py::bytearray();
        PyByteArray_Resize(output.ptr(), output_size);
        assert(size_t(PyByteArray_Size(output.ptr())) == output_size);
        auto output_data = reinterpret_cast<uint8_t *>(PyByteArray_AsString(output.ptr()));

        for (size_t block = 0; block < buf.size / input_block_size; ++block) {
            auto output_block = &output_data[block * output_block_size];
            std::copy_n(&buf.data[block * input_block_size], input_block_size, output_block);
            encode(output_block, input_block_size, &output_block[input_block_size]);
        }

        if (input_remainder > 0) {
            auto output_block = &output_data[output_size - input_remainder - this->ecc_len];
            std::copy_n(&buf.data[buf.size - input_remainder], input_remainder, output_block);
            encode(output_block, input_remainder, &output_block[input_remainder]);
        }

        return output;
    }

    inline bool py_decode(buffer_rw<uint8_t> buf) {
        size_t msg_size = buf.size - ecc_len;
        return decode(buf.data, msg_size, &buf.data[msg_size]);
    }

    inline py::bytearray py_synds(buffer_ro<uint8_t> buf) {
        if (buf.size < ecc_len)
            return {};

        size_t msg_size = buf.size - ecc_len;
        synds_array_t synds_arr;

        synds(buf.data, msg_size, &buf.data[msg_size], synds_arr);

        return py::bytearray(reinterpret_cast<const char *>(synds_arr), ecc_len);
    }

    static inline void register_class(py::module &m) {
        using namespace pybind11::literals;

        py::class_<PyRS256>(m, "RS256")
            .def_property_readonly("ecc_len", [](PyRS256& self) { return self.ecc_len; })

            .def_property("block_len",
                [](PyRS256& self) { return self.default_block_size; },
                &PyRS256::set_default_block_size)

            .def_property_readonly("message_len",
                [](PyRS256& self) { return self.default_block_size - self.ecc_len; })

            .def_property_readonly("gf", [](PyRS256& self) -> auto const& { return self.gf; })

            .def_property_readonly("generator", [](PyRS256& self) {
                return py::bytes(reinterpret_cast<const char *>(self.generator), self.ecc_len + 1); })

            .def_property_readonly("generator_roots", [](PyRS256& self) {
                return py::bytes(reinterpret_cast<const char *>(self.generator_roots), self.ecc_len); })

            .def(py::init<std::optional<uint8_t>, std::optional<uint8_t>, std::optional<uint8_t>, uint8_t, uint16_t>(), R"(
                Instantiate a Reed-Solomon encoder with the given configuration
                )",
                "block_len"_a = py::none(), "message_len"_a = py::none(), "ecc_len"_a = py::none(),
                "primitive"_a = 2, "polynomial"_a = 0x11d)

            .def("__sizeof__", [](PyRS256& self) { return sizeof(self); })

            .def("encode", cast_args(&PyRS256::py_encode),
                R"(Systematic encode)",
                "buffer"_a)

            .def("encode_blocks", cast_args(&PyRS256::py_encode_blocks), R"(
                Encode the input buffer in blocks, storing the parity bytes right next to their corresponding block
                )",
                "buffer"_a, "block_len"_a = py::none())

            .def("decode", cast_args(&PyRS256::py_decode),
                R"(Systematic decode)",
                "buffer"_a)

            .def("_synds", cast_args(&PyRS256::py_synds),
                R"(Compute syndromes)",
                "buffer"_a)

            .doc() = R"(Reed-Solomon coding over :math:`GF(2^8)`)";
    }

private:
    inline uint8_t _get_ecc_len(
            std::optional<uint8_t> block_len,
            std::optional<uint8_t> message_len,
            std::optional<uint8_t> ecc_len) {
        if (ecc_len) {
            return *ecc_len;
        } else if (message_len && block_len) {
            if (*message_len >= *block_len) {
                throw py::value_error("block_len must be greater than message_len");
            }
            return block_len.value_or(default_block_size) - *message_len;
        } else {
            throw py::value_error("Must specify either (block_len, message_len) or ecc_len");
        }
    }
};

