
#  rsi16md.py
#
#  Copyright 2026 Gabriel Machado
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

import libffrs

import ffrs.reference as ref
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted


GF = ref.GF(65537, 1, 3)
ntt = lambda w, x: to_int_list(ref.ntt.ntt(GF, w, rbo_sorted(x)))
intt = lambda w, x: rbo_sorted(to_int_list(ref.ntt.intt(GF, w, x)))


class RSi16md(libffrs.RSi16md):
    def _mix_ecc(self, ecc):
        i = rbo(self.block_len, self.block_len - self.ecc_len)
        w_i = self.gf.pow(self.gf.inv(self.root), i)

        ecc_mix = [self.gf.mul(s, self.gf.sub(0, self.gf.pow(w_i, j))) for j, s in enumerate(ecc)]

        # ecc_root = GF(self.root) ** (self.block_size // self.ecc_size)
        ecc_root = GF(self.gf.exp(self.gf.div(self.gf.log(1), self.ecc_len)))
        ecc_mix = intt(ecc_root, ecc_mix)

        return ecc_mix

    def find_errors(self, msg1: bytearray):
        w = GF(self.root)

        msg1_list = to_int_list(msg1, 2)
        synds = ntt(w, msg1_list)[:self.ecc_len]

        if all(s == 0 for s in synds):
            return {}

        # print("synds:    ", synds)
        synds = GF(synds)

        for err_count in range(self.ecc_len // 2, 0, -1):
            mat = [synds[i:i+err_count] for i in range(err_count)]

            err_loc_coefs = ref.linalg.gaussian_elim(mat, [GF(0) - s for s in synds[err_count:2*err_count]])
            lm = ref.P(GF, [1] + err_loc_coefs[::-1])

            err_pos_rbo = [i for i in range(self.block_len) if int(lm.eval(w.inv().pow(i))) == 0]
            if err_pos_rbo:
                break
        else:
            raise RuntimeError("Decoding failed")

        err_ws = [[w.pow(err_pos_rbo[i] * j) for i in range(err_count)] for j in range(err_count)]
        err_pos = [rbo(self.block_len, i) for i in err_pos_rbo]

        err_mag = ref.linalg.gaussian_elim(err_ws, GF(synds[:err_count]))
        errors = dict(zip(err_pos, map(int, err_mag)))
        print("errors:", errors)
        return errors

    def decode(self, msg1: bytearray):
        errors = self.find_errors(msg1)
        msg1_list = to_int_list(msg1, 2)

        for pos, mag in errors.items():
            msg1_list[pos] = self.gf.sub(msg1_list[pos], mag)

        msg1[:] = to_bytearray(msg1_list, 2)

        return msg1[:self.message_size]
