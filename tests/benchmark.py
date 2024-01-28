
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

import random
import timeit

import ffrs


random.seed(42)

def benchmark(method, ecc_len, block_size, duration=2, update_interval=0.5):
    print(f'\n{method} ecc_len: {ecc_len} block_size: {block_size}')
    RS256 = ffrs.RS256(ecc_len=ecc_len)

    total_time = 0
    total_size = 0
    size = 1 * 1000 * 1000

    while total_time < duration:
        data = bytearray(size)

        time = timeit.timeit(f'RS256.{method}(data, {block_size})',
            number=1, globals=dict(RS256=RS256, data=data))

        throughput = size / time

        print(f'Encode speed: {throughput * 1e-6:.3f} MB/s')

        total_time += time
        total_size += size

        size = round(throughput * update_interval)

    average_throughput = total_size / total_time
    print(f'Average: {average_throughput * 1e-6:.3f} MB/s')


if __name__ == '__main__':
    # benchmark('encode_blocks', 2, 255 - 2)
    # benchmark('encode_blocks', 4, 255 - 4)
    # benchmark('encode_blocks', 6, 255 - 6)
    # benchmark('encode_blocks', 8, 32 - 8)
    benchmark('encode_blocks', 8, 255 - 8)
    # benchmark('encode_blocks', 32, 255 - 32)
    # benchmark('encode_blocks', 64, 255 - 64)
