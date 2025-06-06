import libffrs

from libffrs import GF256, RS256
from libffrs import GFi16, RSi16

import math

import ffrs.reference as ref
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted


GF = ref.GF(65537, 1, 3)
ntt = lambda w, x: to_int_list(ref.ntt.ntt(GF, w, rbo_sorted(x)))
intt = lambda w, x: rbo_sorted(to_int_list(ref.ntt.intt(GF, w, x)))


class RSi16(libffrs.RSi16):
    def encode(self, buf):
        assert 2**round(math.log2(self.ecc_len // 2)) == self.ecc_len // 2, "ecc_len must be a power of 2"
        w = GF(self.roots_of_unity[round(math.log2(self.block_len // 2))])

        # enc_ref = ntt(to_int_list(buf, 2))
        # assert intt(enc_ref) == to_int_list(buf, 2)

        res = libffrs.RSi16.encode(self, buf)

        ecc = to_int_list(res[-self.ecc_len:], 2)
        # print("ws", [int(w.pow(rbo(self.block_len // 2, (self.block_len - self.ecc_len)//2) * j)) for j, s in enumerate(ecc)])
        # print("ecc:     ", ecc)
        # assert enc_ref[:self.ecc_len//2] == ecc

        ecc_mix = self._mix_ecc(w, ecc)

        ecc_unmix = self._unmix_ecc(w, ecc_mix)
        assert ecc_unmix == ecc, "unmix_ecc failed"

        res[-self.ecc_len:] = to_bytearray(ecc_mix, 2)

        return res

    def _mix_ecc(self, w, ecc):
        ecc_mix = [self.gf.div(s, int(w.pow(rbo(self.block_len // 2, (self.block_len - self.ecc_len)//2) * j))) for j, s in enumerate(ecc)]
        # print("ecc_mul: ", ecc_mix)
        ecc_mix = ecc_mix * (self.block_len // self.ecc_len)
        ecc_mix = intt(w, ecc_mix)
        # print("ecc_mix: ", ecc_mix)

        ecc_mix = ecc_mix[:self.ecc_len // 2]
        return ecc_mix

    def _unmix_ecc(self, w, ecc):
        msg_len = (self.block_len - self.ecc_len) // 2
        # print("recv: ", ecc)
        # ecc = [self.gf.mul(s, int(w.pow(rbo(self.block_len//2, (self.block_len - self.ecc_len)//2) * j))) for j, s in enumerate(ecc)]

        ecc_unmix = ecc  + [0] * msg_len
        # print("ecc_unmix: ", ecc_unmix)
        ecc_unmix = ntt(w, ecc_unmix)

        ecc_unmix = ecc_unmix[:self.ecc_len // 2]
        ecc_unmix = [self.gf.mul(s, int(w.pow(rbo(self.block_len // 2, (self.block_len - self.ecc_len)//2) * j))) for j, s in enumerate(ecc_unmix)]
        # print()
        # print("ws", [int(w.pow(rbo(self.block_len // 2, (self.block_len - self.ecc_len)//2) * j)) for j, s in enumerate(ecc)])
        # print("ecc rec: ", ecc_unmix)

        return ecc_unmix

    def find_errors(self, w, synd, n, errors):
        mat = [synd[i:i+errors] for i in range(errors)]

        err_loc_coefs = ref.linalg.gaussian_elim(mat, [GF(0) - s for s in synd[errors:2*errors]])
        lm = ref.P(GF, [1] + err_loc_coefs[::-1])

        roots = [i for i in range(n) if int(lm.eval(w.inv().pow(i))) == 0]
        return roots

    def decode(self, msg1: bytearray):
        size_u16 = self.block_len // 2
        ecc_u16 = self.ecc_len // 2
        w = GF(self.roots_of_unity[round(math.log2(size_u16))])

        msg2 = self.encode(msg1[:-self.ecc_len] + bytearray(self.ecc_len))

        msg1_list = to_int_list(msg1, 2)
        msg2_list = to_int_list(msg2, 2)

        ecc1 = self._unmix_ecc(w, msg1_list[-ecc_u16:])
        ecc2 = self._unmix_ecc(w, msg2_list[-ecc_u16:])

        synds = [self.gf.sub(e2, e1) for e1, e2 in zip(ecc1, ecc2)]
        # print()
        print("synds:    ", synds)
        if all(s == 0 for s in synds):
            return msg1[:self.message_len]

        for err_count in range(1, ecc_u16 // 2 + 1):
            err_pos_rbo = self.find_errors(w, GF(synds), size_u16, err_count)
            if err_pos_rbo:
                break
        else:
            raise RuntimeError("Decoding failed")

        err_ws = [[w.pow(err_pos_rbo[i] * j) for i in range(err_count)] for j in range(err_count)]
        err_pos = [rbo(size_u16, i) for i in err_pos_rbo]

        err_mag = ref.linalg.gaussian_elim(err_ws, GF(synds[:err_count]))
        print(size_u16, "errors:", dict(zip(err_pos, map(int, err_mag))))

        for pos, mag in zip(err_pos, map(int, err_mag)):
            msg1_list[pos] = self.gf.sub(msg1_list[pos], mag)

        # TODO: if pos > msg_len: remix ecc

        msg1[:] = to_bytearray(msg1_list, 2)

        return msg1[:self.message_len]
