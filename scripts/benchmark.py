#  benchmark.py
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

import argparse
import csv
import platform
import random
import re
import subprocess
import time
import timeit

import ffrs

import matplotlib.pyplot as plt


def get_system_info():
    cpu = None

    if platform.system() == "Windows":
        command = "wmic cpu get name"
        output = subprocess.check_output(command, shell=True).decode()
        cpu = output.split("\n")[1].strip()
    elif platform.system() == "Linux":
        command = r"grep '^\s*model name' /proc/cpuinfo"
        output = subprocess.check_output(command, shell=True).decode()
        for line in output.split("\n"):
            cpu = re.sub(r".*model name.*:\s*(.*)", r"\1", line, 1).strip()
            break

    info = platform.platform()
    if cpu:
        info = "\n".join([info, cpu])
    return info


def random_bytearray(size):
    chunk_size = 2**20
    if size <= chunk_size:
        data = bytearray(random.randbytes(size))
    else:
        data = bytearray(size)
        for i in range(0, size, chunk_size):
            data[i : i + chunk_size] = random.randbytes(min(chunk_size, size - i))
        remaining = size % chunk_size
        if remaining:
            data[-remaining:] = random.randbytes(remaining)
    return data


def benchmark_throughput(
    statement,
    env,
    input_size,
    setup=None,
    cpu_cache_flush_size=50 * 2**20,  # 50 MB
    duration=5,
    update_interval=1,
    cooldown=5,
):
    number = 1
    repeat = 5

    assert "random" not in env
    env["random"] = random

    best_time = float("inf")

    # some cool down time seems to make the test results more repeatable
    while cooldown > 0:
        print(f"Cooldown {cooldown:2}", end="\r")
        cooldown -= 1
        time.sleep(1)

    elapsed_time = 0
    while elapsed_time < duration:

        start_time = time.time()

        times = timeit.repeat(
            statement,
            # operation on a large chunk of memory to hopefully flush CPU cache
            setup=f"random.randbytes({cpu_cache_flush_size}); {setup if setup else ""}",
            number=number,
            repeat=repeat,
            globals=env,
        )

        best_time = min(best_time, min(times))
        throughput = input_size * number / min(times)
        print(f"Speed: {throughput * 1e-6:.3f} MB/s (repeat: {repeat})")

        end_time = time.time()
        elapsed_time += end_time - start_time

        repeat = max(1, round(repeat * update_interval / (end_time - start_time)))

    peak_throughput = input_size / best_time
    print(f"Peak: {peak_throughput * 1e-6:.3f} MB/s")
    return peak_throughput


def benchmark_ber(rs, tries=100):
    if rs.message_size <= 500e6:
        data = random.randbytes(rs.message_size)
    else:
        data = bytearray(rs.message_size)
        chunk_size = 2**20
        for i in range(0, rs.message_size, chunk_size):
            data[i : i + chunk_size] = random.randbytes(min(chunk_size, rs.message_size - i))
        remaining = rs.message_size % chunk_size
        if remaining:
            data[-remaining:] = random.randbytes(remaining)

    ecc = rs.encode(data)
    decode_result = []

    for error_rate in range(50, 200, 25):
        error_count = round(error_rate / 10000 * rs.block_size)
        successes = 0
        for j in range(tries):
            error_positions = random.sample(range(rs.block_size), error_count)
            # print(f"\rerror_rate: {error_rate} ({error_count}/{rs.block_size}) attempt: {j:3}/{tries}", end="")

            data_t = bytearray(data)
            ecc_t = bytearray(ecc)

            for i in error_positions:
                if i < rs.message_size:
                    data_t[i] ^= random.randint(1, 255)
                else:
                    ecc_t[i - rs.message_size] ^= random.randint(1, 255)

            rs.repair(data_t, ecc_t)

            successes += int(data == data_t and ecc == ecc_t)
        decode_result.append((error_count, successes, tries))

    print("\r\033[2K", end="")
    for error_count, successes, tries in decode_result:
        success_rate = successes / tries
        print(f"Error rate: {error_count / rs.block_size:.2} {round(error_count)}, Success rate: {success_rate:.1%}")
    print()

    return decode_result


def enc256_benchmark(method, ecc_len, block_size):
    print(f"\n{method}block_size: {block_size} ecc_len: {ecc_len} ")
    RS256 = ffrs.RS256(block_size, ecc_len=ecc_len)

    return (f"rs.{method}(data)", dict(rs=RS256))


def run_enc_benchmarks(config):
    table = []

    try:
        with open("benchmark.csv", mode="r") as file:
            reader = csv.reader(file)
            for row in reader:
                table.append(row)
    except FileNotFoundError:
        pass

    if table:
        assert len(table) == len(config) + 1
    else:
        table = [["method", "ecc_len", "block_size"]]
        for cfg in config:
            table.append(list(map(str, cfg)))

    header = table[0]

    throughput = []
    for method, ecc_len, block_size in config:
        thr = benchmark_throughput(enc256_benchmark(method, ecc_len, block_size))
        throughput.append(f"{thr * 1e-6:.2f} MB/s")

    if ffrs.compiler_info in header:
        index = header.index(ffrs.compiler_info)
        for row, cfg, thr in zip(table[1:], config, throughput):
            assert row[:3] == list(map(str, cfg))
            row[index] = thr
    else:
        header.append(ffrs.compiler_info)
        for row, cfg, thr in zip(table[1:], config, throughput):
            assert row[:3] == list(map(str, cfg))
            row.append(thr)

    with open("benchmark.csv", mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(table)


def parse_size(size_str: str) -> int:
    size_str = size_str.strip().lower()

    multipliers = {
        "k": 1024,
        "m": 1024**2,
        "g": 1024**3,
        "t": 1024**4,
    }

    if size_str[-1] in multipliers:
        number = float(size_str[:-1]) if "." in size_str[:-1] else int(size_str[:-1])
        return int(number * multipliers[size_str[-1]])
    else:
        return int(size_str)


def parse_ratio(ratio_str: str) -> tuple[float, float]:
    ratio_str = ratio_str.strip().lower()

    if "/" in ratio_str:
        parts = ratio_str.split("/")
    elif ":" in ratio_str:
        parts = ratio_str.split(":")
    else:
        raise ValueError(f"Invalid ratio format: {ratio_str}")

    if len(parts) != 2:
        raise ValueError(f"Invalid ratio format: {ratio_str}")

    numerator = parse_size(parts[0])
    denominator = parse_size(parts[1])

    return numerator, denominator


def circ(args):
    block_size = parse_size(args.block_size)

    ecc_ratio_num, ecc_ratio_den = parse_ratio(args.ecc_ratio)
    assert ecc_ratio_num == 1
    assert ecc_ratio_den > 0 and (ecc_ratio_den & (ecc_ratio_den - 1)) == 0

    hash_ratio_num, hash_ratio_den = parse_ratio(args.hash_ratio)

    block_size_mul = block_size
    block_size_power = 0
    while block_size_mul > 0 and block_size_mul & 1 == 0:
        block_size_mul >>= 1
        block_size_power += 1

    print(
        f"block_size: {block_size_mul} * 2**{block_size_power} ecc_ratio: 1/{ecc_ratio_den} hash_ratio: {hash_ratio_num}/{hash_ratio_den}"
    )

    ecc_ratio_num *= 2
    ecc_ratio_den *= 2

    results = []

    # TODO: validate RSi16md with 64k block
    while block_size % (hash_ratio_den * ecc_ratio_den) == 0 and ecc_ratio_den <= 32768:
        outer_interleave = block_size // (hash_ratio_den * ecc_ratio_den)
        rs = ffrs.CIRC(hash_ratio_den // 2, hash_ratio_num // 2, ecc_ratio_den, ecc_ratio_num, outer_interleave)
        print(rs)
        assert rs.block_size == block_size, rs.block_size

        # benchmark_throughput((f"rs.encode(data)", dict(rs=rs)), input_size=rs.message_size)
        res = benchmark_ber(rs)
        results.append((ecc_ratio_den, outer_interleave, res))

        ecc_ratio_num *= 2
        ecc_ratio_den *= 2

    cmap = plt.cm.Blues
    for i, (rso_block_len, outer_interleave, res) in enumerate(results):
        success_rates = [successes / tries for _, successes, tries in res]
        error_counts = [error_count for error_count, _, _ in res]
        color = cmap((i + 1) / (len(results) + 1))
        plt.plot(
            error_counts,
            success_rates,
            color=color,
            label=f"rso_block_size: {rso_block_len * 2} outer_interleave: {outer_interleave}",
        )

    plt.ylabel("repair success rate")
    plt.xlabel("errors")
    plt.xscale("log", base=2)
    plt.legend()
    plt.tight_layout()
    plt.show()


re_algo = re.compile(r"(RSi16md|CIRC)\((\w+=)?\d+(,\s*(\w+=)?\d+)*\)")


def parse_algo(algo):
    assert re_algo.match(algo), f"Invalid format: {algo}"
    return eval("ffrs." + algo)


def cmd_throughput(args):
    rs = parse_algo(args.algo)
    print(get_system_info())
    print(rs.message_size / 2**20, rs.ecc_size / 2**20, rs.ecc_size / rs.message_size)
    data = random_bytearray(rs.message_size)
    benchmark_throughput(f"rs.encode(data)", dict(rs=rs, data=data), input_size=rs.message_size)


def decode_throughput(args):
    block_size = parse_size(args.block_size)

    ecc_ratio_num, ecc_ratio_den = parse_ratio(args.ecc_ratio)
    assert ecc_ratio_num == 1
    assert ecc_ratio_den > 0 and (ecc_ratio_den & (ecc_ratio_den - 1)) == 0

    extra_args = dict()
    if args.no_simd:
        extra_args["simd_x4"] = False
        extra_args["simd_x8"] = False
        extra_args["simd_x16"] = False

    rsi_ecc, rsi_block = parse_ratio(args.hash_ratio)
    rso_block = block_size // rsi_block
    rs = ffrs.CIRC(
        rsi_block // 2,
        rsi_ecc // 2,
        rso_block,
        rso_block * ecc_ratio_num // ecc_ratio_den,
        args.outer_interleave,
        **extra_args,
    )

    assert rs.block_size == block_size * args.outer_interleave

    print(get_system_info())
    print(rs)
    print(rs.message_size / 2**20, rs.ecc_size / 2**20, rs.ecc_size / rs.message_size)

    data = random_bytearray(rs.message_size)

    if args.errors_per_col is None:
        args.errors_per_col = rs.outer_ecc_len // 2

    def add_errors(data):
        for col in range(rs.inner_message_len):
            error_positions = random.sample(range(rs.outer_message_len), args.errors_per_col)
            for pos in error_positions:
                data[pos * rs.inner_message_size + col * 2] ^= random.randint(1, 255)
                data[pos * rs.inner_message_size + col * 2 + 1] ^= random.randint(1, 255)

    ecc = rs.encode(data)
    # orig = bytearray(data)

    benchmark_throughput(
        f"rs.repair(data, ecc)",
        dict(rs=rs, data=data, ecc=ecc, add_errors=add_errors),
        setup="add_errors(data)",
        input_size=len(data),
    )

    # assert data == orig


def main():
    parser = argparse.ArgumentParser(description="Benchmark")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Throughput
    parser_thr = subparsers.add_parser("throughput", help="Benchmark throughput")
    parser_thr.add_argument("algo")
    parser_thr.set_defaults(func=cmd_throughput)

    # # BER
    # parser_ber = subparsers.add_parser("ber", help="Remove an item")
    # parser_ber.add_argument("algo")
    # parser_ber.set_defaults(func=cmd_ber)

    # CIRC
    parser_circ = subparsers.add_parser("circ", help="Benchmark CIRC decoding throughput")
    parser_circ.add_argument("block_size")
    parser_circ.add_argument("ecc_ratio")
    parser_circ.add_argument("hash_ratio")
    parser_circ.set_defaults(func=circ)

    # CIRC decode throughput
    parser_dec = subparsers.add_parser("dec_throughput", help="Benchmark CIRC decode throughput")
    parser_dec.add_argument("--no-simd", action="store_true", help="Disable SIMD optimizations")
    parser_dec.add_argument("--errors-per-col", type=int, default=None, help="Number of errors to introduce per column")
    parser_dec.add_argument("block_size")
    parser_dec.add_argument("ecc_ratio")
    parser_dec.add_argument("hash_ratio")
    parser_dec.add_argument("outer_interleave", nargs="?", default=1, type=int)
    parser_dec.set_defaults(func=decode_throughput)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
