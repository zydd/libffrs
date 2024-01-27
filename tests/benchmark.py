
#  test_lib_rs256.py
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

RS256 = ffrs.RS256(ecc_len=255 - 223)

random.seed(42)

def test_encode_blocks():
    total_time = 0
    for i in range(0, 1000, 10):
        size = i * 1000 * 1000
        data = bytearray(size)

        time = timeit.timeit("RS256.encode_blocks(data, 223)",
            number=1, globals=dict(RS256=RS256, data=data))

        throughput = size * 1e-6 / time

        print(f'Encode speed: {throughput:.3f} MB/s')

        total_time += time
        if total_time > 5:
            break

if __name__ == '__main__':
    test_encode_blocks()
