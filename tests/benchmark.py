
#  benchmark.py
#
#  Copyright 2024 Gabriel Machado
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
        command = 'wmic cpu get name'
        output = subprocess.check_output(command, shell=True).decode()
        cpu = output.split('\n')[1].strip()
    elif platform.system() == "Linux":
        command = r"grep '^\s*model name' /proc/cpuinfo"
        output = subprocess.check_output(command, shell=True).decode()
        for line in output.split('\n'):
            cpu = re.sub(r'.*model name.*:\s*(.*)', r'\1', line, 1).strip()
            break

    info = platform.platform()
    if cpu:
        info = '\n'.join([info, cpu])
    return info


def benchmark_peak(
        method, ecc_len, block_size,
        input_size=100 * 1000 * 1000,  # 100 MB
        cpu_cache_flush_size=50 * 1000 * 1000,  # 50 MB
        duration=5, update_interval=1, cooldown=5):
    print(f'\n{method} ecc_len: {ecc_len} block_size: {block_size}')
    RS256 = ffrs.RS256(ecc_len=ecc_len)

    number = 1
    repeat = 5

    best_time = float('inf')

    # some cool down time seems to make the test results more repeatable
    time.sleep(cooldown)

    data = random.randbytes(input_size)

    elapsed_time = 0
    while elapsed_time < duration:
        start_time = time.time()

        times = timeit.repeat(
            f'rs.{method}(data, {block_size})',

            # operation in a large block of memory to (hopefully) flush CPU cache
            setup=f'random.randbytes({cpu_cache_flush_size})',

            number=number, repeat=repeat,
            globals=dict(rs=RS256, data=data, random=random))

        best_time = min(best_time, min(times))
        throughput = input_size * number / min(times)
        print(f'Encode speed: {throughput * 1e-6:.3f} MB/s (repeat: {repeat})')

        end_time = time.time()
        elapsed_time += end_time - start_time

        repeat = max(1, round(repeat * update_interval / (end_time - start_time)))

    peak_throughput = input_size / best_time
    print(f'Peak: {peak_throughput * 1e-6:.3f} MB/s')
    return peak_throughput


if __name__ == '__main__':
    print(get_system_info())

    config = [
        # ('method', ecc_len, block_size)

        # Max block size
        ('encode_blocks', 2, 255 - 2),
        ('encode_blocks', 4, 255 - 4),
        ('encode_blocks', 6, 255 - 6),
        ('encode_blocks', 8, 255 - 8),
        ('encode_blocks', 10, 255 - 10),
        ('encode_blocks', 12, 255 - 12),
        ('encode_blocks', 16, 255 - 16),
        ('encode_blocks', 24, 255 - 24),
        ('encode_blocks', 32, 255 - 32),
        ('encode_blocks', 64, 255 - 64),
        ('encode_blocks', 128, 255 - 128),

        # CD-ROM CIRC
        ('encode_blocks', 32 - 28, 28),  # inner
        ('encode_blocks', 32 - 24, 24),  # outer

        # DVD CIRC
        ('encode_blocks', 208 - 192, 192),  # inner
        ('encode_blocks', 182 - 172, 172),  # outer
    ]

    throughput = []
    for method, ecc_len, block_size in config:
        if len(sys.argv) == 1 or str(ecc_len) in sys.argv[1:]:
            throughput.append(benchmark_peak(method, ecc_len, block_size))
