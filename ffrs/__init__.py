import libffrs

from libffrs import GF256, RS256
from libffrs import GFi16, RSi16

import functools
import math

import ffrs.reference as ref
from ffrs.reference.util import to_int_list, to_bytearray, rbo

GF = ref.GF(65537, 1, 3)


class RSi16(libffrs.RSi16):
    def find_errors(self, w, synd, n, errors):
        mat = [synd[i:i+errors] for i in range(errors)]

        err_loc_coefs = ref.linalg.gaussian_elim(mat, [GF(0) - s for s in synd[errors:2*errors]])
        lm = ref.P(GF, [1] + err_loc_coefs[::-1])

        roots = [i for i in range(n) if int(lm.eval(w.inv().pow(i))) == 0]
        return roots

    def decode(self, msg1: bytearray):
        size_u16 = self.block_len // 2
        ecc_u16 = self.ecc_len // 2

        msg2 = msg1[:-self.ecc_len] + bytearray(self.ecc_len)
        self.encode(msg2)

        msg1_list = to_int_list(msg1, 2)
        msg2_list = to_int_list(msg2, 2)

        ecc1 = msg1_list[-ecc_u16:]
        ecc2 = msg2_list[-ecc_u16:]

        w = GF(self.roots_of_unity[round(math.log2(size_u16))])
        synds = [self.gf.sub(e2, e1) for e1, e2 in zip(ecc1, ecc2)]
        if all(s == 0 for s in synds):
            return

        for err_count in range(1, ecc_u16 // 2 + 1):
            err_pos_rbo = self.find_errors(w, GF(synds), size_u16, err_count)
            if err_pos_rbo:
                break
        else:
            raise RuntimeError("Decoding failed")

        err_ws = [[w.pow(err_pos_rbo[i] * j) for i in range(err_count)] for j in range(err_count)]
        err_pos = [rbo(size_u16, i) for i in err_pos_rbo]

        err_mag = ref.linalg.gaussian_elim(err_ws, GF(synds[:err_count]))

        for pos, mag in zip(err_pos, map(int, err_mag)):
            msg1_list[pos] = self.gf.sub(msg1_list[pos], mag)

        msg1[:] = to_bytearray(msg1_list, 2)

        return pos
