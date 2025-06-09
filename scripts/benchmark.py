
#  benchmark.py
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

import csv
import platform
import random
import re
import subprocess
import sys
import time
import timeit

import ffrs

random.seed(42)


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


def benchmark_throughput(
        benchmark_config,
        input_size=100 * 2**20,  # 100 MB
        cpu_cache_flush_size=50 * 2**20,  # 50 MB
        duration=5, update_interval=1, cooldown=5
    ):

    number = 1
    repeat = 5

    best_time = float("inf")

    # some cool down time seems to make the test results more repeatable
    while cooldown > 0:
        print(f"Cooldown {cooldown:2}", end="\r")
        cooldown -= 1
        time.sleep(1)

    data = random.randbytes(input_size)

    elapsed_time = 0
    while elapsed_time < duration:
        code, env = benchmark_config
        start_time = time.time()

        env.update(dict(
            data=data,
            random=random,
        ))

        times = timeit.repeat(
            code,
            # operation on a large chunk of memory to hopefully flush CPU cache
            setup=f"random.randbytes({cpu_cache_flush_size})",
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


def ntt_benchmark(block_size):
    print(f"\nntt block_size: {block_size * 2} bytes")
    GF = ffrs.GFi16(65537, 3)
    w = next(i for i in range(2, GF.prime) if i**block_size % GF.prime == 1)

    return (f"gf.ntt_blocks16(data, {w}, {block_size})", dict(gf=GF))

def show_help():
    print("benchmark.py [enc|ntt]")


def main():
    print(get_system_info())

    config = [
        # ("method", ecc_len, block_size)

        # Max block size
        ("encode_blocks", 2, 255),
        ("encode_blocks", 4, 255),
        ("encode_blocks", 6, 255),
        ("encode_blocks", 8, 255),
        ("encode_blocks", 10, 255),
        ("encode_blocks", 12, 255),
        ("encode_blocks", 14, 255),
        ("encode_blocks", 16, 255),
        ("encode_blocks", 20, 255),
        ("encode_blocks", 24, 255),
        ("encode_blocks", 32, 255),
        ("encode_blocks", 64, 255),
        ("encode_blocks", 128, 255),
    ]

    if len(sys.argv) < 2:
        return show_help()

    block_size = 1024
    ecc_len = block_size // 8

    fn = sys.argv[1]

    if fn == "enc":
        if len(sys.argv) > 2:
            for arg in sys.argv[2:]:
                ecc_len = int(arg)
                benchmark_throughput(enc256_benchmark("encode_blocks", ecc_len, 255))
        else:
            run_enc_benchmarks(config)
    elif fn == "ntt":
        benchmark_throughput(ntt_benchmark(1024//2), input_size=1 * 2**20)
    elif fn == "enci16":
        benchmark_throughput((f"rs.encode_blocks(data)", dict(rs=ffrs.RSi16(block_size, ecc_len=ecc_len))), input_size=(block_size - ecc_len) * 2**10)
    elif fn == "enci16v":
        benchmark_throughput((f"rs.encode_blocks(data)", dict(rs=ffrs.RSi16md(block_size, ecc_len=ecc_len))), input_size=(block_size - ecc_len) * 2**10)
    else:
        return show_help()


if __name__ == "__main__":
    main()
