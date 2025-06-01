
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
import functools

import ffrs
import ffrs.reference
import ffrs.reference.ntt


from common import to_int_list, to_bytearray

# random.seed(42)

gfref = ffrs.reference.GF(65537, 1, 3)
GF = ffrs.GFi32(65537, 3)

def rbo(max, i):
    nbits = round(math.log2(max))
    return int(bin(i)[2:].zfill(nbits)[::-1], 2)

def rbo_sorted(a):
    assert len(a) and len(a) & (len(a) - 1) == 0, "len(a) must be a power of 2"
    # return list(a)
    return [a[rbo(len(a), i)] for i in range(len(a))]


@pytest.mark.parametrize('rs', [
    ffrs.RSi16(8, ecc_len=2*2),
    ffrs.RSi16(32, ecc_len=4*2),
    ffrs.RSi16(256, ecc_len=8*2),
    ffrs.RSi16(256, ecc_len=9*2),
    ffrs.RSi16(256, ecc_len=10*2),
    ffrs.RSi16(512, ecc_len=4*2),
])
class TestRS:
    def test_encode(self, rs):
        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        buf = [random.randrange(0, 2**16) for _ in range(size_u16 - ecc_u16)] + [0] * ecc_u16

        w = gfref(rs.roots_of_unity[round(math.log2(size_u16))])

        ref = to_int_list(ffrs.reference.ntt.ntt(gfref, gfref(w), rbo_sorted(buf)))
        res_i = to_int_list(ffrs.reference.ntt.intt(gfref, gfref(w), ref))

        res = to_bytearray(buf, 2)
        rs.encode(res)
        res = to_int_list(res, 2)

        res2 = rs.gf.ntt16(to_bytearray(buf, 2), w)
        res2 = to_int_list(res2, 2)

        assert rbo_sorted(res_i) == buf
        assert ref[:ecc_u16] == res[-ecc_u16:]

    def test_encode_decode(self, rs):
        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        orig = [random.randrange(0, 2**16) for _ in range(size_u16 - ecc_u16)] + [0] * ecc_u16

        msg_enc = to_bytearray(orig, 2)
        rs.encode(msg_enc)

        msg_enc_err = bytearray(msg_enc)
        assert msg_enc == msg_enc_err
        msg_enc_err = to_int_list(msg_enc_err, 2)

        err_vec = [0] * size_u16

        error_positions = dict()
        while len(error_positions) < ecc_u16//2:
            # Do not corrupt codeword for now
            i = random.randrange(size_u16 - ecc_u16)
            if i in error_positions:
                continue
            error_positions[i] = random.randrange(2**16)
            msg_enc_err[i] = rs.gf.add(msg_enc_err[i], error_positions[i])
            err_vec[i] = rs.gf.add(err_vec[i], error_positions[i])

        print()
        print("errors:", error_positions)

        assert msg_enc_err != msg_enc

        msg_enc_err = to_bytearray(msg_enc_err, 2)
        res_pos = _find_err_pos(rs, msg_enc_err, orig, err_vec, err_count=len(error_positions))
        assert set(res_pos) == set(error_positions.keys())

        # TODO: correct errors

    def test_encode_blocks(self, rs):
        blocks = 3

        size_u16 = rs.block_len // 2
        ecc_u16 = rs.ecc_len // 2
        msg_u16 = rs.message_len // 2
        buf = [random.randrange(1, 2**16) for _ in range(msg_u16 * blocks)]

        w = gfref(rs.roots_of_unity[round(math.log2(size_u16))])

        buf_ntt = [to_int_list(ffrs.reference.ntt.ntt(gfref, gfref(w), rbo_sorted(buf[msg_u16 * blk:msg_u16 * (blk + 1)] + [0] * ecc_u16))) for blk in range(blocks)]
        res_i = [to_int_list(ffrs.reference.ntt.intt(gfref, gfref(w), rbo_sorted(buf_ntt[blk]))) for blk in range(blocks)]

        res = rs.encode_blocks(to_bytearray(buf, 2))
        res = to_int_list(res, 2)

        for i in range(blocks):
            start = (i + 1) * size_u16 - ecc_u16
            assert buf_ntt[i][:ecc_u16] == res[start:start + ecc_u16]


from ffrs.reference import P
from ffrs.reference.linalg import gaussian_elim

def find_errors(w, synd, n, errors):
    mat = [synd[i:i+errors] for i in range(errors)]

    err_loc_coefs = gaussian_elim(mat, [gfref(0) - s for s in synd[errors:2*errors]])
    lm = P(gfref, [1] + err_loc_coefs[::-1])

    roots = [i for i in range(n) if int(lm.eval(w.inv().pow(i))) == 0]
    return roots

def _find_err_pos(rs, msg1, orig, err_vec, err_count):
    size_u16 = rs.block_len // 2
    ecc_u16 = rs.ecc_len // 2

    msg2 = msg1[:-rs.ecc_len] + bytearray(rs.ecc_len)
    rs.encode(msg2)

    msg1 = to_int_list(msg1, 2)
    msg2 = to_int_list(msg2, 2)

    ecc1 = msg1[-ecc_u16:]
    ecc2 = msg2[-ecc_u16:]

    w = gfref(rs.roots_of_unity[round(math.log2(size_u16))])
    synds_ref = [functools.reduce(type(w).__add__, 
                                  [w.pow(j*i) * gfref(e) for i, e in enumerate(rbo_sorted(err_vec))])
                                    for j in range(ecc_u16)]


    synds = [rs.gf.sub(e2, e1) for e1, e2 in zip(ecc1, ecc2)]
    print("synds   :", synds)
    assert synds == synds_ref

    # breakpoint()

    pos = find_errors(w, gfref(synds), size_u16, err_count)
    pos = [rbo(size_u16, i) for i in pos]
    print("pos     :", pos)

    return pos
