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

from .backup import maximize_interleaving
from . import cli
from .cli import log


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


def main(args):
    input_size = args.input_size.get()
    rs = maximize_interleaving(args, input_size)
    log.debug("platform: %s", platform.platform())
    log.debug("processor: %s", platform.processor())
    log.debug(
        "python: %s %s %s %s", platform.python_implementation(), platform.python_version(), *platform.architecture()
    )
    log.debug("ffrs compiled with: %s", ffrs.compiler_info)
    log.info("codec: %s", rs)

    match args.function.get():
        case "encode":
            log.info("buffer size: %s", ffrs.util.format_size(rs.message_size))
            data = ffrs.create_buffer(rs.message_size)
            benchmark_throughput(f"rs.encode(data)", dict(rs=rs, data=data), input_size=rs.message_size)

        case other:
            raise ffrs.par.FfrsParException(f"unknown function: {other}")

    log.info("done")
    return 0


def arg_parser(parser):
    cli_parser = cli.CLI(parser, main)
    cli_parser.parser.add_argument("function", metavar="func")
    cli_parser.parser.add_argument(
        "--input-size",
        metavar="size",
        type=cli.parse_size,
        default=cli.DEFAULT("100m", cli.parse_size),
        help="Input size for benchmarking",
    )
    return cli_parser


if __name__ == "__main__":
    cli.main()
