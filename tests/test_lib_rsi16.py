
#  test_lib_rsi16.py
#
#  Copyright 2025 Gabriel Machado
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import math
import pytest
import random

import ffrs
import ffrs.reference
import ffrs.reference.ntt
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted

# random.seed(42)

ref_gf = ffrs.reference.GF(65537, 1, 3)

ref_ntt_ct = lambda w, buf: to_int_list(ffrs.reference.ntt.ntt(ref_gf, ref_gf(w), rbo_sorted(buf)))
ref_intt_ct = lambda w, buf: rbo_sorted(to_int_list(ffrs.reference.ntt.intt(ref_gf, ref_gf(w), buf)))
ref_ntt_gs = lambda w, buf: rbo_sorted(to_int_list(ffrs.reference.ntt.ntt(ref_gf, ref_gf(w), buf)))
ref_intt_gs = lambda w, buf: to_int_list(ffrs.reference.ntt.intt(ref_gf, ref_gf(w), rbo_sorted(buf)))

ref_ntt, ref_intt = ref_ntt_ct, ref_intt_ct
# ref_ntt, ref_intt = ref_ntt_gs, ref_intt_gs


@pytest.mark.parametrize('rs', [
    # ffrs.RSi16(8, ecc_len=2*2),
    ffrs.RSi16(32, ecc_len=4*2),
    ffrs.RSi16(256, ecc_len=4*2),
    ffrs.RSi16(256, ecc_len=8*2),
    ffrs.RSi16(256, ecc_len=16*2),
    ffrs.RSi16(512, ecc_len=4*2),
])
class TestRS:
    def test_encode(self, rs):
        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        buf = [random.randrange(0, 2**16) for _ in range(size_u16 - ecc_u16)] + [0] * ecc_u16

        w = ref_gf(rs.roots_of_unity[round(math.log2(size_u16))])

        ref = ref_ntt(w, buf)
        res_i = ref_intt(w, ref)

        res = rs.encode(to_bytearray(buf, 2))
        res = to_int_list(res, 2)

        ecc_mix = ref[:ecc_u16]
        ecc_mix = [rs.gf.sub(0, s) for s in ecc_mix]
        ecc_mix = [rs.gf.div(s, int(w.pow(rbo(rs.block_len // 2, (rs.block_len - rs.ecc_len)//2) * j))) for j, s in enumerate(ecc_mix)]
        ecc_mix = ecc_mix * (size_u16 // ecc_u16)
        # print("ecc_mix:", ecc_mix)
        ecc_mix = ref_intt(w, ecc_mix)
        ref_ecc = ecc_mix[:ecc_u16]
        # print("ref_ecc:", ref_ecc)

        assert res_i == buf
        assert ref_ecc == res[-ecc_u16:]

    def test_encode_decode(self, rs):
        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        orig = [random.randrange(0, 2**16) for _ in range(size_u16 - ecc_u16)] + [0] * ecc_u16

        msg_enc = to_int_list(rs.encode(to_bytearray(orig, 2)), 2)

        msg_enc_err = list(msg_enc)

        assert rs.find_errors(to_bytearray(msg_enc, 2)) == {}


        error_positions = dict()

        # Always add one error to the codeword
        if len(error_positions) < max(ecc_u16//2, 1):
            i = random.randrange(size_u16 - ecc_u16, size_u16)
            error_positions[i] = random.randrange(2**16)
            msg_enc_err[i] = rs.gf.add(msg_enc_err[i], error_positions[i])

        while len(error_positions) < ecc_u16//2:
            i = random.randrange(size_u16)
            if i in error_positions:
                continue
            error_positions[i] = random.randrange(2**16)
            msg_enc_err[i] = rs.gf.add(msg_enc_err[i], error_positions[i])

        w = ref_gf(rs.roots_of_unity[round(math.log2(size_u16))])
        print()
        print("errs: ", error_positions)

        assert msg_enc_err != msg_enc

        msg_enc_err_dec = rs.decode(to_bytearray(msg_enc_err, 2))

        assert orig[:-ecc_u16] == to_int_list(msg_enc_err_dec, 2)
        assert error_positions == rs.find_errors(to_bytearray(msg_enc_err, 2))

    def test_encode_blocks(self, rs):
        blocks = 3

        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        msg_u16 = rs.message_len // 2
        buf = [random.randrange(1, 2**16) for _ in range(msg_u16 * blocks)]

        w = ref_gf(rs.roots_of_unity[round(math.log2(size_u16))])

        buf_ntt = [ref_ntt(w, buf[msg_u16 * blk:msg_u16 * (blk + 1)] + [0] * ecc_u16) for blk in range(blocks)]

        res = rs.encode_blocks(to_bytearray(buf, 2))
        res = to_int_list(res, 2)

        for i in range(blocks):
            start = (i + 1) * size_u16 - ecc_u16
            assert buf_ntt[i][:ecc_u16] == res[start:start + ecc_u16]

