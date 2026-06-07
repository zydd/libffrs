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
import ffrs.par
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


def main(args):
    log.debug("platform: %s", platform.platform())
    log.debug("processor: %s", platform.processor())
    log.debug(
        "python: %s %s %s %s", platform.python_implementation(), platform.python_version(), *platform.architecture()
    )
    log.debug("ffrs compiled with: %s", ffrs.compiler_info)

    match args.function.get():
        case "encode":
            benchmark_encode(args)

        case other:
            raise ffrs.par.FfrsParException(f"unknown function: {other}")

    log.info("done")
    return 0


def arg_parser(parser):
    # Add module group first
    benchmark = parser.add_argument_group("benchmark")

    cli_parser = cli.CLI(parser, main)
    cli_parser.parser.add_argument("function", metavar="func")

    benchmark.add_argument(
        "-s",
        "--input-size",
        metavar="size",
        type=cli.parse_size,
        default=cli.DEFAULT("100m", cli.parse_size),
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
    return cli_parser


if __name__ == "__main__":
    cli.main()
