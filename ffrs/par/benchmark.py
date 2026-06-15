#  par/benchmark.py
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

"""Benchmark function"""

import platform
import random
import time
import timeit

import ffrs
import ffrs.util

from . import cli
from .opt import maximize_interleaving
from .cli import log


def benchmark_throughput(
    args,
    statement,
    env,
    input_size,
    setup=None,
):
    duration = args.duration.get()
    cpu_cache_flush_size = args.cpu_cache_flush_size.get()
    update_interval = args.update_interval.get()
    cooldown = args.cooldown.get()
    number = 1
    repeat = args.min_repeats.get()

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


def fill_random(buffer):
    chunk_size = 2**20  # 2 MiB
    if len(buffer) <= chunk_size:
        buffer[:] = random.randbytes(len(buffer))
    else:
        for i in range(0, len(buffer), chunk_size):
            buffer[i : i + chunk_size] = random.randbytes(min(chunk_size, len(buffer) - i))
        remaining = len(buffer) % chunk_size
        if remaining:
            buffer[-remaining:] = random.randbytes(remaining)


def benchmark_encode(args):
    input_size = args.input_size.get()
    rs = maximize_interleaving(args, input_size)
    log.info("codec: %s", rs)
    log.info("input size: %s", ffrs.util.format_size(input_size))
    log.info(
        "buffer size: %s (+%.1f%%)", ffrs.util.format_size(rs.message_size), (rs.message_size / input_size - 1) * 100
    )
    data = ffrs.create_buffer(rs.message_size)
    benchmark_throughput(args, f"rs.encode(data)", dict(rs=rs, data=data), input_size=rs.message_size)


def benchmark_decode(args):
    input_size = args.input_size.get()
    rs = maximize_interleaving(args, input_size)
    log.info("codec: %s", rs)
    log.info("input size: %s", ffrs.util.format_size(input_size))
    log.info(
        "buffer size: %s (+%.1f%%)", ffrs.util.format_size(rs.message_size), (rs.message_size / input_size - 1) * 100
    )
    data = ffrs.create_buffer(rs.message_size)
    fill_random(data)
    ecc = rs.encode(data)

    errors_per_col = rs.outer_ecc_len // 2

    def add_errors(data):
        for col in range(rs.inner_message_len):
            error_positions = random.sample(range(rs.outer_block_len), errors_per_col)
            for pos in error_positions:
                if pos < rs.outer_message_len:
                    data[pos * rs.inner_message_size + col * 2] ^= random.randint(1, 255)
                    data[pos * rs.inner_message_size + col * 2 + 1] ^= random.randint(1, 255)
                else:
                    pos -= rs.outer_message_len
                    ecc[pos * rs.inner_message_size + col * 2] ^= random.randint(1, 255)
                    ecc[pos * rs.inner_message_size + col * 2 + 1] ^= random.randint(1, 255)

    benchmark_throughput(
        args,
        f"rs.repair(data, ecc)",
        dict(rs=rs, data=data, ecc=ecc, add_errors=add_errors),
        setup="add_errors(data)",
        input_size=rs.message_size,
    )


def benchmark_ber(args):
    input_size = args.input_size.get()
    tries = args.ber_tries.get()
    rs = maximize_interleaving(args, input_size)
    log.info("input size: %s", ffrs.util.format_size(input_size))
    log.info("tries: %s", tries)

    data = ffrs.create_buffer(rs.message_size)
    fill_random(data)

    results = []
    while rs.outer_ecc_len >= 16:
        log.info("codec: %s", rs)

        assert rs.message_size == len(data)
        ecc = rs.encode(data)

        decode_result = []

        no_success = False
        for error_rate in range(50, 10000, 25):
            error_count = round(error_rate * rs.block_size / 10000)
            log.debug("testing error rate: %.3f (errors: %s)", error_rate / 10000, error_count)
            successes = 0
            t0 = time.time()
            for j in range(tries):
                error_positions = random.sample(range(rs.block_size), error_count)

                data_t = bytearray(data)
                ecc_t = bytearray(ecc)

                for i in error_positions:
                    if i < rs.message_size:
                        data_t[i] ^= random.randint(1, 255)
                    else:
                        ecc_t[i - rs.message_size] ^= random.randint(1, 255)

                rs.repair(data_t, ecc_t)

                successes += int(data == data_t and ecc == ecc_t)
            t1 = time.time()

            success_rate = successes / tries
            decode_result.append((error_count, successes, tries))
            throughput = input_size * tries / (t1 - t0)
            log.info(
                "input error rate: %.3f (%d), repair success rate: %.1f%% throughput: %s/s",
                error_count / rs.block_size,
                round(error_count),
                success_rate * 100,
                ffrs.util.format_size(throughput),
            )

            if successes == 0:
                if no_success:
                    break
                no_success = True

        results.append((rs, decode_result))
        rs = ffrs.CIRC16(
            rs.inner_block_len,
            rs.inner_ecc_len,
            rs.outer_block_len // 2,
            rs.outer_ecc_len // 2,
            rs.interleave * 2,
            primitive=rs.primitive,
            simd_x4=rs.simd_x4,
            simd_x8=rs.simd_x8,
            simd_x16=rs.simd_x16,
        )

    if args.show_plot.get():
        import matplotlib.pyplot as plt
        import matplotlib.ticker

        cmap = plt.cm.Blues
        for i, (rs, res) in enumerate(results):
            success_rates = [successes / tries for _, successes, tries in res]
            error_rates = [error_count / rs.block_size for error_count, _, _ in res]
            color = cmap((len(results) - i + 1) / (len(results) + 1))
            plt.plot(
                error_rates,
                success_rates,
                color=color,
                label=str(rs),
            )

        plt.ylabel("repair success rate")
        plt.xlabel("error rate")
        plt.gca().xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1))
        # plt.xscale("log", base=2)
        plt.legend()
        plt.tight_layout()
        plt.show()


FUNCTIONS = {
    "decode": benchmark_decode,
    "encode": benchmark_encode,
    "ber": benchmark_ber,
}


def main(args):
    log.debug("platform: %s", platform.platform())
    log.debug("processor: %s", platform.processor())
    log.debug(
        "python: %s %s %s %s", platform.python_implementation(), platform.python_version(), *platform.architecture()
    )
    log.debug("ffrs compiled with: %s", ffrs.compiler_info)

    FUNCTIONS[args.function.get()](args)

    log.info("done")
    return 0


def create_parser(parent=None):
    parser = cli.create_parser(__loader__.name, __doc__, parent)
    from .common_args import add_common_args

    add_common_args(parser)

    benchmark = parser.parser.add_argument_group("benchmark")

    benchmark.add_argument("function", metavar="func", choices=FUNCTIONS.keys(), help="Function to benchmark")
    benchmark.add_argument(
        "-s",
        "--input-size",
        metavar="size",
        type=cli.parse_size,
        default=cli.DEFAULT("128m", cli.parse_size),
        help="Input size for benchmarking",
    )
    benchmark.add_argument(
        "-t",
        "--duration",
        metavar="time",
        type=int,
        default=cli.DEFAULT(5),
        help="Benchmark duration",
    )
    benchmark.add_argument(
        "-c",
        "--cooldown",
        metavar="time",
        type=int,
        default=cli.DEFAULT(5),
        help="Cooldown before starting benchmark",
    )
    benchmark.add_argument(
        "--min-repeats",
        metavar="n",
        type=int,
        default=cli.DEFAULT(5),
        help="Minimum number of executions",
    )
    benchmark.add_argument(
        "--cpu-cache-flush-size",
        metavar="size",
        type=cli.parse_size,
        default=cli.DEFAULT("50m", cli.parse_size),
        help="Operation on a large chunk of memory to mitigate CPU caching effects on measured performance",
    )
    benchmark.add_argument(
        "--update-interval",
        metavar="time",
        type=int,
        default=cli.DEFAULT(1),
        help="Interval to print current throughput to stdout",
    )
    benchmark.add_argument(
        "--ber-tries",
        metavar="n",
        type=int,
        default=cli.DEFAULT(100),
        help="Number of tries for bit error rate benchmark",
    )
    benchmark.add_argument(
        "--show-plot",
        action="store_true",
        help="Show plot",
    )
    return parser


if __name__ == "__main__":
    cli.main()
