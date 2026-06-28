#  test_lib_rsi16.py
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

import pytest
import random

import ffrs
import ffrs.reference as ref
import ffrs.reference.rs as ref_rs
import ffrs.reference.ntt
from ffrs.reference.util import to_int_list, to_bytearray, rbo, rbo_sorted, randbytes

ref_gf = ref.GF(65537, 1, 3)

ref_ntt = lambda w, buf: to_int_list(ffrs.reference.ntt.ntt(ref_gf, ref_gf(w), rbo_sorted(buf)))
ref_intt = lambda w, buf: rbo_sorted(to_int_list(ffrs.reference.ntt.intt(ref_gf, ref_gf(w), buf)))


def add_aligned_errors(rs: ffrs.RSi16, buf, ecc, n):
    rows = list(random.sample(range(rs.rs_block_len), n))
    for row in rows:
        if row < rs.rs_message_len:
            for col in range(rs.interleave):
                # message
                offset = 2 * rs.message_offset(row, col)
                buf[offset] ^= random.randint(1, 255)
                buf[offset + 1] ^= random.randint(1, 255)
        else:
            for col in range(rs.interleave):
                # ecc
                offset = 2 * rs.ecc_offset(row - rs.rs_message_len, col)
                ecc[offset] ^= random.randint(1, 255)
                ecc[offset + 1] ^= random.randint(1, 255)
    return rows


def add_unaligned_errors(rs: ffrs.RSi16, msg, ecc, n):
    error_positions = []
    for col in range(rs.interleave):
        error_positions.append([])
        for row in random.sample(range(rs.rs_block_len), n):
            error_positions[-1].append(row)
            if row < rs.rs_message_len:
                # message
                offset = 2 * rs.message_offset(row, col)
                msg[offset] ^= random.randint(1, 255)
                msg[offset + 1] ^= random.randint(1, 255)
            else:
                # ecc
                offset = 2 * rs.ecc_offset(row - rs.rs_message_len, col)
                ecc[offset] ^= random.randint(1, 255)
                ecc[offset + 1] ^= random.randint(1, 255)
    return error_positions


class BaseTestRS:
    def test_encode(self, rs: ffrs.RSi16):
        buf = randbytes(rs.message_size)

        res = rs.encode(buf)

        w = ref_gf(rs.root)
        buf_ntt = ref_ntt(w, to_int_list(buf + bytearray(rs.ntt_size - len(buf))))

        ecc_mix = ref_rs.mix_ecc(rs, ref_gf, buf_ntt[: rs.ecc_len])

        assert to_bytearray(ecc_mix) == res

    @pytest.mark.parametrize("grace", [0, 1, 2, 3, 4])
    def test_repair_unknown_locations(self, rs: ffrs.RSi16, grace):
        msg_orig = randbytes(rs.message_size)
        ecc_orig = rs.encode(msg_orig)
        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        _error_positions = add_aligned_errors(rs, msg_err, ecc_err, max(rs.ecc_len // 2 - grace, 1))

        assert msg_err != msg_orig or ecc_err != ecc_orig

        res = rs.repair(msg_err, ecc_err)
        assert res == ffrs.RepairStatus.RepairOk or all(s == 0 for s in rs._synd(msg_err, ecc_err))  # FIXME

        assert msg_err == msg_orig
        assert ecc_err == ecc_orig

    @pytest.mark.parametrize("grace", [0, 1, 2, 3, 4])
    def test_repair_known(self, rs: ffrs.RSi16, grace):
        msg_orig = randbytes(rs.message_size)
        ecc_orig = rs.encode(msg_orig)
        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        error_positions = add_aligned_errors(rs, msg_err, ecc_err, max(rs.ecc_len - grace, 1))

        assert msg_err != msg_orig or ecc_err != ecc_orig

        res = rs.repair(msg_err, ecc_err, error_positions)
        assert res == ffrs.RepairStatus.RepairOk

        assert msg_err == msg_orig
        assert ecc_err == ecc_orig

    def test_repair_no_errors(self, rs: ffrs.RSi16):
        msg_orig = randbytes(rs.message_size)
        ecc_orig = rs.encode(msg_orig)
        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        res = rs.repair(msg_err, ecc_err, [])
        assert res == ffrs.RepairStatus.NoErrors

        assert msg_err == msg_orig
        assert ecc_err == ecc_orig

    def test_repair_unknown_no_errors(self, rs: ffrs.RSi16):
        msg_orig = randbytes(rs.message_size)
        ecc_orig = rs.encode(msg_orig)
        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        res = rs.repair(msg_err, ecc_err)
        assert res == ffrs.RepairStatus.RepairOk or all(s == 0 for s in rs._synd(msg_err, ecc_err))  # FIXME

        assert msg_err == msg_orig
        assert ecc_err == ecc_orig

    def test_encode_blocks_empty(self, rs: ffrs.RSi16):
        assert bytearray() == rs.encode(bytearray())

    def test_encode_blocks_single(self, rs: ffrs.RSi16):
        buf = randbytes(rs.message_size)

        buf_enc = rs.encode(buf)
        buf_enc_blk = rs.encode(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_size
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize("count", [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 100])
    def test_encode_blocks_multiple(self, rs: ffrs.RSi16, count):
        buf = randbytes(rs.message_size * count)

        buf_enc = [rs.encode(buf[i * rs.message_size : (i + 1) * rs.message_size]) for i in range(count)]
        buf_enc = b"".join(buf_enc)
        buf_enc_blk = rs.encode(buf)

        assert len(buf_enc) == len(buf_enc_blk) == rs.ecc_size * count
        assert buf_enc == buf_enc_blk

    @pytest.mark.parametrize("interleave", list(range(2, 15)) + [32, 48, 100, 256])
    def test_encode_interleaved(self, rs: ffrs.RSi16, interleave):
        assert rs.interleave == 1
        rsi = ffrs.RSi16(
            rs.block_len,
            rs.ecc_len,
            interleave=interleave,
            simd_x4=rs.simd_x4,
            simd_x8=rs.simd_x8,
            simd_x16=rs.simd_x16,
        )

        data = list(range(interleave * rs.message_len))
        buf = to_bytearray(data)

        interleaved_chunk = rsi.encode(buf)

        for i in range(interleave):
            buf_enc = rs.encode(to_bytearray(data[i::interleave]))

            # Sequential ecc
            # assert interleaved_chunk[i * rs.ecc_size:(i + 1) * rs.ecc_size] == buf_enc

            # Interleaved ecc
            assert to_bytearray(to_int_list(interleaved_chunk)[i::interleave]) == buf_enc

    @pytest.mark.parametrize("interleave", list(range(1, 16)) + [32, 48, 50])
    @pytest.mark.parametrize("grace", range(4))
    def test_repair_interleaved_unknown_unaligned(self, rs: ffrs.RSi16, interleave, grace):
        assert rs.interleave == 1
        rsi = ffrs.RSi16(
            rs.block_len,
            rs.ecc_len,
            interleave=interleave,
            simd_x4=rs.simd_x4,
            simd_x8=rs.simd_x8,
            simd_x16=rs.simd_x16,
        )

        msg_orig = randbytes(rsi.message_size)
        ecc_orig = rsi.encode(msg_orig)

        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        _error_positions = add_unaligned_errors(rs, msg_err, ecc_err, max(rs.ecc_len // 2 - grace, 1))

        msg_err0 = to_bytearray(to_int_list(msg_err)[::interleave])
        ecc_err0 = to_bytearray(to_int_list(ecc_err)[::interleave])
        res = rs.repair(msg_err0, ecc_err0)
        assert res == ffrs.RepairStatus.RepairOk or all(s == 0 for s in rs._synd(msg_err0, ecc_err0))  # FIXME

        assert msg_err != msg_orig or ecc_err != ecc_orig
        res = rsi.repair(msg_err, ecc_err)
        assert msg_err == msg_orig
        assert ecc_err == ecc_orig
        assert res == ffrs.RepairStatus.RepairOk

    @pytest.mark.parametrize("interleave", list(range(1, 16)) + [32, 33, 48, 50])
    @pytest.mark.parametrize("grace", range(4))
    def test_repair_interleaved_unknown_aligned(self, rs: ffrs.RSi16, interleave, grace):
        assert rs.interleave == 1
        rsi = ffrs.RSi16(
            rs.block_len,
            rs.ecc_len,
            interleave=interleave,
            simd_x4=rs.simd_x4,
            simd_x8=rs.simd_x8,
            simd_x16=rs.simd_x16,
        )

        msg_orig = randbytes(rsi.message_size)
        ecc_orig = rsi.encode(msg_orig)

        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        add_aligned_errors(rsi, msg_err, ecc_err, max(rsi.rs_ecc_len // 2 - grace, 1))

        assert msg_err != msg_orig or ecc_err != ecc_orig
        rsi.repair(msg_err, ecc_err)
        assert msg_err == msg_orig
        assert ecc_err == ecc_orig

    @pytest.mark.parametrize("interleave", list(range(1, 16)) + [32, 33, 48, 50])
    @pytest.mark.parametrize("grace", range(4))
    def test_repair_interleaved_known_aligned(self, rs: ffrs.RSi16, interleave, grace):
        assert rs.interleave == 1
        rsi = ffrs.RSi16(
            rs.block_len,
            rs.ecc_len,
            interleave=interleave,
            simd_x4=rs.simd_x4,
            simd_x8=rs.simd_x8,
            simd_x16=rs.simd_x16,
        )

        msg_orig = randbytes(rsi.message_size)
        ecc_orig = rsi.encode(msg_orig)

        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        rows = add_aligned_errors(rsi, msg_err, ecc_err, max(rsi.rs_ecc_len - grace, 1))

        assert msg_err != msg_orig or ecc_err != ecc_orig
        rsi.repair(msg_err, ecc_err, rows)
        assert msg_err == msg_orig
        assert ecc_err == ecc_orig

    def test_sugiyama_no_errors(self, rs: ffrs.RSi16, subtests):
        locator, evaluator = rs._sugiyama(bytearray(rs.ecc_size))
        assert locator[0] == 0

    @pytest.mark.parametrize("n", [4, 8, 16])
    def test_sugiyama_no_errors_n(self, rs: ffrs.RSi16, n, subtests):
        sugiyama = getattr(rs, f"_sugiyama{n}")

        locators, evaluators = sugiyama(bytearray(n * rs.ecc_size))

        for i in range(n):
            locator = locators[i * rs.ecc_len : (i + 1) * rs.ecc_len]
            assert locator[0] == 0

    def test_sugiyama(self, rs: ffrs.RSi16, subtests):
        for err_count in range(1, rs.ecc_len // 2 + 1):
            with subtests.test(err_count=err_count):
                buf = bytearray(rs.message_size)
                ecc = bytearray(rs.ecc_size)

                error_locations = add_aligned_errors(rs, buf, ecc, err_count)

                locator_ref = ref_rs.locator(ref_gf, rs.root, map(rs.ntt.rbo, error_locations))
                locator_ref = list(map(int, locator_ref.x)) + [0] * (rs.ecc_len - len(locator_ref.x))

                synd = rs._synd(buf, ecc)
                locator, evaluator = rs._sugiyama(synd)
                roots = rs._roots(locator)

                roots_ref = rs._roots(locator_ref)

                assert set(error_locations) == set(roots)
                assert roots == roots_ref
                assert locator == locator_ref

    @pytest.mark.parametrize("n", [4, 8, 16])
    def test_sugiyama_aligned_n(self, rs: ffrs.RSi16, n, subtests):
        sugiyama = getattr(rs, f"_sugiyama{n}")

        for err_count in range(1, rs.ecc_len // 2 + 1):
            with subtests.test(err_count=err_count):
                buf = bytearray(rs.message_size * n)
                ecc = bytearray(rs.ecc_size * n)

                error_locations = random.sample(range(rs.message_len), err_count)
                for i in range(n):
                    for e in error_locations:
                        buf[i * rs.message_size + e * 2] = random.randrange(1, 256)

                synd = rs._synd(buf, ecc)
                locators, evaluators = sugiyama(synd)

                for i in range(n):
                    start = i * rs.ecc_len
                    locator_ref, evaluator_ref = rs._sugiyama(synd[start : start + rs.ecc_len])
                    locator = locators[start : start + rs.ecc_len]
                    evaluator = evaluators[start : start + rs.ecc_len]
                    assert set(error_locations) == set(rs._roots(locator))
                    assert locator == locator_ref
                    assert evaluator == evaluator_ref

    @pytest.mark.parametrize("n", [4, 8, 16])
    def test_sugiyama_max_unaligned_n(self, rs: ffrs.RSi16, n, subtests):
        sugiyama = getattr(rs, f"_sugiyama{n}")

        for err_count in range(1, rs.ecc_len // 2 + 1):
            with subtests.test(err_count=err_count):
                buf = bytearray(rs.message_size * n)
                errors = []
                for i in range(n):
                    errors.append(random.sample(range(rs.message_len), err_count))
                    for e in errors[-1]:
                        buf[i * rs.message_size + e * 2] = random.randrange(1, 256)

                synd = rs._synd(buf, bytearray(n * rs.ecc_size))
                locators, evaluators = sugiyama(synd)

                for i in range(n):
                    start = i * rs.ecc_len
                    locator_ref, evaluator_ref = rs._sugiyama(synd[start : start + rs.ecc_len])
                    locator = locators[start : start + rs.ecc_len]
                    evaluator = evaluators[start : start + rs.ecc_len]
                    assert locator == locator_ref
                    assert evaluator == evaluator_ref
                    assert set(errors[i]) == set(rs._roots(locator))

    @pytest.mark.parametrize("n", [4])
    def test_sugiyama_unaligned_random_n(self, rs: ffrs.RSi16, n, subtests):
        sugiyama = getattr(rs, f"_sugiyama{n}")

        for err_count in range(1, rs.ecc_len // 2 + 1):
            with subtests.test(err_count=err_count):
                buf = bytearray(rs.message_size * n)
                errors = []
                for i in range(n):
                    errors.append(random.sample(range(rs.message_len), random.randrange(1, err_count + 1)))
                    for e in errors[-1]:
                        buf[i * rs.message_size + e * 2] = random.randrange(1, 256)

                synd = rs._synd(buf, bytearray(n * rs.ecc_size))

                locator_refs = []
                evaluator_refs = []
                for i in range(n):
                    start = i * rs.ecc_len
                    # capsys.readouterr()
                    locator_ref, evaluator_ref = rs._sugiyama(synd[start : start + rs.ecc_len])
                    # with open(f"log{i}.txt", "w") as f:
                    #     f.write(capsys.readouterr().out)
                    locator_refs.append(locator_ref)
                    evaluator_refs.append(evaluator_ref)

                # capsys.readouterr()
                locators, evaluators = sugiyama(synd)
                # with open(f"log{4}.txt", "w") as f:
                #     f.write(capsys.readouterr().out)

                for i in range(n):
                    start = i * rs.ecc_len
                    locator = locators[start : start + rs.ecc_len]
                    evaluator = evaluators[start : start + rs.ecc_len]
                    assert locator == locator_refs[i]
                    assert evaluator == evaluator_refs[i]
                    assert set(errors[i]) == set(rs._roots(locator))


@pytest.mark.parametrize(
    "rs",
    [
        ffrs.RSi16(4, ecc_len=2),
        ffrs.RSi16(8, ecc_len=2),
        ffrs.RSi16(16, ecc_len=4),
        ffrs.RSi16(16, ecc_len=8),
        ffrs.RSi16(128, ecc_len=2),
        ffrs.RSi16(256, ecc_len=32),
        # ffrs.RSi16(1024, ecc_len=128),
        # ffrs.RSi16(4096, ecc_len=512),
        ffrs.RSi16(4, ecc_len=2, simd_x16=False),
        ffrs.RSi16(16, ecc_len=4, simd_x16=False),
        ffrs.RSi16(128, ecc_len=2, simd_x16=False),
        ffrs.RSi16(256, ecc_len=32, simd_x16=False),
        # ffrs.RSi16(1024, ecc_len=128, simd_x16=False),
        # ffrs.RSi16(4096, ecc_len=512, simd_x16=False),
        ffrs.RSi16(4, ecc_len=2, simd_x16=False, simd_x8=False),
        ffrs.RSi16(16, ecc_len=4, simd_x16=False, simd_x8=False),
        ffrs.RSi16(128, ecc_len=2, simd_x16=False, simd_x8=False),
        ffrs.RSi16(256, ecc_len=32, simd_x16=False, simd_x8=False),
        # ffrs.RSi16(1024, ecc_len=128, simd_x16=False, simd_x8=False),
        # ffrs.RSi16(4096, ecc_len=512, simd_x16=False, simd_x8=False),
        ffrs.RSi16(4, ecc_len=2, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.RSi16(16, ecc_len=4, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.RSi16(128, ecc_len=2, simd_x16=False, simd_x8=False, simd_x4=False),
        ffrs.RSi16(256, ecc_len=32, simd_x16=False, simd_x8=False, simd_x4=False),
        # ffrs.RSi16(1024, ecc_len=128, simd_x16=False, simd_x8=False, simd_x4=False),
        # ffrs.RSi16(4096, ecc_len=512, simd_x16=False, simd_x8=False, simd_x4=False),
    ],
)
class TestRSPower2(BaseTestRS):
    pass


@pytest.mark.parametrize(
    "rs",
    [
        ffrs.RSi16(16 * 13, ecc_len=8),
        ffrs.RSi16(256 * 3, ecc_len=8),
    ],
)
class TestRSMult(BaseTestRS):
    pass


@pytest.mark.parametrize(
    "rs",
    [
        ffrs.RSi16(4096, ecc_len=8),
        ffrs.RSi16(4096, ecc_len=256),
        ffrs.RSi16(2048, ecc_len=128),
        ffrs.RSi16(4096, ecc_len=2048),
        ffrs.RSi16(32768, ecc_len=256),
        ffrs.RSi16(65536, ecc_len=256),
    ],
)
class TestRSLarge:
    def test_repair_no_errors(self, rs: ffrs.RSi16):
        msg = randbytes(rs.message_size)
        ecc = rs.encode(msg)
        msg_err = bytearray(msg)
        ecc_err = bytearray(ecc)

        res = rs.repair(msg_err, ecc_err, [])
        assert res == ffrs.RepairStatus.NoErrors

        assert msg_err == msg
        assert ecc_err == ecc

    def test_repair_unknown_no_errors(self, rs: ffrs.RSi16):
        msg = randbytes(rs.message_size)
        ecc = rs.encode(msg)
        msg_err = bytearray(msg)
        ecc_err = bytearray(ecc)

        res = rs.repair(msg_err, ecc_err)
        assert res == ffrs.RepairStatus.RepairOk or all(s == 0 for s in rs._synd(msg_err, ecc_err))  # FIXME

        assert msg_err == msg
        assert ecc_err == ecc

    @pytest.mark.parametrize("grace", range(4))
    def test_repair_unknown_locations(self, rs: ffrs.RSi16, grace):
        msg_orig = randbytes(rs.message_size)
        ecc_orig = rs.encode(msg_orig)

        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        _error_positions = add_aligned_errors(rs, msg_err, ecc_err, max(rs.ecc_len // 2 - grace, 1))

        assert msg_err != msg_orig or ecc_err != ecc_orig

        res = rs.repair(msg_err, ecc_err)
        assert res == ffrs.RepairStatus.RepairOk or all(s == 0 for s in rs._synd(msg_err, ecc_err))  # FIXME

        assert msg_err == msg_orig
        assert ecc_err == ecc_orig

    @pytest.mark.parametrize("grace", range(4))
    def test_repair_known_locations(self, rs: ffrs.RSi16, grace):
        msg_orig = randbytes(rs.message_size)
        ecc_orig = rs.encode(msg_orig)

        msg_err = bytearray(msg_orig)
        ecc_err = bytearray(ecc_orig)

        error_positions = add_aligned_errors(rs, msg_err, ecc_err, max(rs.ecc_len - grace, 1))

        assert msg_err != msg_orig or ecc_err != ecc_orig

        res = rs.repair(msg_err, ecc_err, error_positions)
        assert res == ffrs.RepairStatus.RepairOk or all(s == 0 for s in rs._synd(msg_err, ecc_err))  # FIXME

        assert msg_err == msg_orig
        assert ecc_err == ecc_orig
