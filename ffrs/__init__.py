import libffrs

from libffrs import *

import math

import ffrs.reference as ref
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted


GF = ref.GF(65537, 1, 3)
ntt = lambda w, x: to_int_list(ref.ntt.ntt(GF, w, rbo_sorted(x)))
intt = lambda w, x: rbo_sorted(to_int_list(ref.ntt.intt(GF, w, x)))


class RSi16(libffrs.RSi16):
    def _mix_ecc(self, w, ecc):
        ecc_mix = ecc
        i = rbo(self.block_size // 2, (self.block_size - self.ecc_len) // 2)
        w_i = self.gf.pow(self.gf.inv(w), i)

        ecc_mix = [self.gf.mul(s, self.gf.pow(w_i, j)) for j, s in enumerate(ecc_mix)]

        # print("ecc_mul: ", ecc_mix)
        ecc_mix = ecc_mix * (self.block_size // self.ecc_len)
        ecc_mix = intt(w, ecc_mix)
        # print("ecc_mix: ", ecc_mix)

        ecc_mix = ecc_mix[:self.ecc_len // 2]
        ecc_mix = [self.gf.sub(0, s) for s in ecc_mix]
        return ecc_mix

    def find_errors(self, msg1: bytearray):
        size_u16 = self.block_size // 2
        ecc_u16 = self.ecc_len // 2
        w = GF(self.roots_of_unity[round(math.log2(size_u16))])

        # msg2 = self.encode(msg1[:-self.ecc_len] + bytearray(self.ecc_len))

        msg1_list = to_int_list(msg1, 2)
        # msg2_list = to_int_list(msg2, 2)

        # ecc1 = self._unmix_ecc(w, msg1_list[-ecc_u16:])
        # ecc2 = self._unmix_ecc(w, msg2_list[-ecc_u16:])

        # synds = [self.gf.sub(e1, e2) for e1, e2 in zip(ecc1, ecc2)]
        # # print()
        # print("synds:    ", synds)
        synds = ntt(w, msg1_list)[:ecc_u16]
        if all(s == 0 for s in synds):
            return {}
        synds = GF(synds)

        for err_count in range(ecc_u16 // 2, 0, -1):
            mat = [synds[i:i+err_count] for i in range(err_count)]

            err_loc_coefs = ref.linalg.gaussian_elim(mat, [GF(0) - s for s in synds[err_count:2*err_count]])
            lm = ref.P(GF, [1] + err_loc_coefs[::-1])

            err_pos_rbo = [i for i in range(size_u16) if int(lm.eval(w.inv().pow(i))) == 0]
            if err_pos_rbo:
                break
        else:
            raise RuntimeError("Decoding failed")

        err_ws = [[w.pow(err_pos_rbo[i] * j) for i in range(err_count)] for j in range(err_count)]
        err_pos = [rbo(size_u16, i) for i in err_pos_rbo]

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

        return msg1[:self.message_len]
