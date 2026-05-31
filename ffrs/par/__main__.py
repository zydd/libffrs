#  par/__main__.py
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

import base64
import hashlib
import logging
import math
import os
import re
import sys
import time

import ffrs
import ffrs.par
import ffrs.util

from . import constraints
from .cli import CLI
from .solver import Solver


class ConfigurationError(ffrs.par.FfrsParException):
    pass


class MaxLevelFilter(logging.Filter):
    def __init__(self, level):
        self.level = level

    def filter(self, record):
        return record.levelno < self.level


class ColorFormatter(ffrs.par.ColorFormatter):
    format_string = "%(levelname)s: "

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.FORMATS[logging.INFO] = logging.Formatter("%(message)s")


logger = logging.getLogger("__main__")
logger.setLevel(logging.INFO)

stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setLevel(logging.DEBUG)
stdout_handler.addFilter(MaxLevelFilter(logging.WARNING))
stdout_handler.setFormatter(ColorFormatter())
logger.addHandler(stdout_handler)

stderr_handler = logging.StreamHandler(sys.stderr)
stderr_handler.setLevel(logging.WARNING)
stderr_handler.setFormatter(ColorFormatter())
logger.addHandler(stderr_handler)


def format_size(n):
    units = ["", "k", "m", "g", "t", "p"]
    colors = [
        "\033[36m",  # cyan
        "\033[32m",  # green
        "\033[33m",  # yellow
        "\033[31m",  # red
        "\033[35m",  # magenta
        "\033[94m",  # light blue
    ]
    if n <= 0:
        return str(n)
    parts = []
    for eng in range(len(units) - 1, -1, -1):
        shift = eng * 10
        value = 1 << shift
        if n >= value:
            q, n = divmod(n, value)
            if q:
                parts.append(f"{colors[eng]}{q}{units[eng]}\033[0m")
    return "".join(parts)


def _ansi_ljust(s, w):
    return s + " " * max(0, w - len(re.sub(r"\033\[[0-9;]*m", "", s)))


def print_config(config):
    logger.debug("config: {")
    for k in sorted(config):
        val = config[k]
        if val is None:
            val = "\033[35m[[[...]]]\033[0m"
        elif len(val) > 32:
            val = f"\033[35m[[[{len(val)}]]]\033[0m"
        else:
            val = "[" + ", ".join((map(format_size, sorted(val)))) + "]"

        logger.debug(f"  {repr(k):23}: {val},")
    logger.debug("}\n")


def print_config_args(config):
    assert (
        len(config["inner_block"])
        == len(config["inner_ecc"])
        == len(config["outer_block"])
        == len(config["outer_ecc"])
        == len(config["interleave"])
        == 1
    )
    logger.debug(
        f"CIRC16({next(iter(config['inner_block'])) // 2}, {next(iter(config['inner_ecc'])) // 2},"
        f" {next(iter(config['outer_block']))}, {next(iter(config['outer_ecc']))}, {next(iter(config['interleave']))})"
    )


def instantiate_config(config):
    return ffrs.CIRC16(
        next(iter(config["inner_block"])) // 2,
        next(iter(config["inner_ecc"])) // 2,
        next(iter(config["outer_block"])),
        next(iter(config["outer_ecc"])),
        next(iter(config["interleave"])),
    )


def branch_config(solver, config, arg, value):
    t0 = time.time()
    # logger.debug("branch %s = %s", arg, value)
    # config2 = copy.deepcopy(config)  # slower than making a new dict
    config2 = {k: set(v) if v is not None else None for k, v in config.items()}
    config2[arg] = [value]
    solver.solve(config2, [arg])
    logger.debug("branching %s = %s took %.2f ms", arg, value, (time.time() - t0) * 1000)
    return config2


def optimize(solver, mode, config, arg):
    t0 = time.time()
    if mode == "max":
        arg_values = sorted(config[arg], reverse=True)
    elif mode == "min":
        arg_values = sorted(config[arg])
    else:
        raise NotImplementedError

    if len(config[arg]) > 1:
        for i, val in enumerate(arg_values):
            try:
                config2 = branch_config(solver, config, arg, val)
                logger.debug("optimizing %s = %s (index: %d) took %.2f ms", arg, val, i, (time.time() - t0) * 1000)
                return config2
            except ValueError:
                pass
        raise ConfigurationError(f"could not optimize {arg}")
    return config


def maximize_interleaving(args, file_size):
    t0 = time.time()
    file_size = max(file_size, 8192)

    config = dict(
        block=args.block.get(),
        message=None,  # choose based on file_size
        ecc=args.ecc.get(),
        inner_block=args.inner_block.get(range(2, 65536 + 1, 2)),
        inner_message=range(2, 65536 + 1, 2),  # choose based on file_size
        inner_ecc=args.inner_ecc.get(set(2**i for i in range(1, 16))),
        outer_block=args.outer_block.get(range(2, 65536 + 1, 2)),
        outer_message=args.outer_message.get(range(2, 65536 + 1, 2)),
        outer_ecc=args.outer_ecc.get(set(2**i for i in range(1, 16))),
        outer_ecc_ratio_num=(
            set(num for num, den in args.ecc_ratio.get()) if args.ecc_ratio.has_value() else range(1, 1 + 1)
        ),
        outer_ecc_ratio_den=(
            set(den for num, den in args.ecc_ratio.get()) if args.ecc_ratio.has_value() else range(1, 1024 + 1)
        ),
        interleave=args.interleave.get(),
        inner_interleaved_ecc=args.outer_interleaved_ecc.get(),
        outer_interleaved_ecc=args.outer_interleaved_ecc.get(),
        outer_interleave=args.outer_interleave.get(),
    )

    free_variables = {k for k in config if k not in args or not getattr(args, k).has_value()}
    free_variables |= {"message", "inner_message"}
    if args.ecc_ratio.has_value():
        free_variables.remove("outer_ecc_ratio_den")
        free_variables.remove("outer_ecc_ratio_num")

    print_config(config)

    constraints.append(f"{file_size} <= {{message}} <= {2 * file_size}")

    solver = Solver(config.keys(), constraints, free_variables)

    try:
        solver.solve(config)
    except ValueError:
        logger.exception("could not determine configuration")
        # return
        raise

    for outer_ecc in config["outer_ecc"]:
        try:
            logger.debug("-" * 72)
            logger.debug("branch outer_ecc")
            config2 = branch_config(solver, config, "outer_ecc", outer_ecc)
            print_config(config2)

            logger.debug("-" * 72)
            logger.debug("compute outer_interleave")
            max_padding = 0.05 if file_size <= 16000 else 0.01
            config2["outer_interleave"] = [
                math.ceil(file_size / outer_message) + i
                for outer_message in config2["outer_message"]
                for i in range(min(100, math.ceil(file_size * max_padding / outer_message)))
            ]
            solver.solve(config2, ["outer_interleave"])
            print_config(config2)

            logger.debug("-" * 72)
            logger.debug("maximize inner_message")
            config2 = optimize(solver, "max", config2, "inner_message")
            print_config(config2)

            logger.debug("-" * 72)
            logger.debug("minimize inner_ecc")
            config2 = optimize(solver, "min", config2, "inner_ecc")
            print_config(config2)

            logger.debug("-" * 72)
            logger.debug("minimize message")
            config2 = optimize(solver, "min", config2, "message")
            print_config(config2)

            print("solver iterations:", solver._constraint_evaluations)
            print("max message ratio:", max(config2["message"]) / file_size)
            print("max hash rato", max(config2["inner_interleaved_ecc"]) / file_size)
            config = config2
            break
        except ValueError:
            logger.debug("could not set outer_ecc to %s", outer_ecc)
    else:
        raise ConfigurationError("could not optimize interleaving")

    print_config(config)
    t1 = time.time() - t0
    if t1 > 10e-3:
        logger.warning("time to solve: %.2f [ms]", t1 * 1000)
    return instantiate_config(config)


def circ(args):
    config = dict(
        block=args.block.get(),
        message=args.message.get(),
        ecc=args.ecc.get(),
        inner_block=args.inner_block.get(range(2, 65536 + 1, 2)),
        inner_message=args.inner_message.get(range(2, 65536 + 1, 2)),
        inner_ecc=args.inner_ecc.get(set(2**i for i in range(1, 16))),
        outer_block=args.outer_block.get(range(2, 65536 + 1, 2)),
        outer_message=args.outer_message.get(range(2, 65536 + 1, 2)),
        outer_ecc=args.outer_ecc.get(set(2**i for i in range(1, 16))),
        outer_ecc_ratio_num=(
            set(num for num, den in args.ecc_ratio.get()) if args.ecc_ratio.has_value() else range(1, 1 + 1)
        ),
        outer_ecc_ratio_den=(
            set(den for num, den in args.ecc_ratio.get()) if args.ecc_ratio.has_value() else range(1, 1024 + 1)
        ),
        interleave=args.interleave.get(),
        inner_interleaved_ecc=args.outer_interleaved_ecc.get(),
        outer_interleaved_ecc=args.outer_interleaved_ecc.get(),
        outer_interleave=args.outer_interleave.get(),
    )

    free_variables = {k for k in config if k not in args or not getattr(args, k).has_value()}
    if args.ecc_ratio.has_value():
        free_variables.remove("outer_ecc_ratio_den")
        free_variables.remove("outer_ecc_ratio_num")

    print_config(config)

    solver = Solver(config.keys(), constraints, free_variables)

    try:
        solver.solve(config)
    except ValueError:
        logger.exception("could not determine configuration")
        return

    print_config(config)

    domain_size = solver.domain_size(config)
    logger.debug("domain size: %s", domain_size)

    def maximize(arg):
        nonlocal domain_size
        if not getattr(args, arg).is_set() and domain_size > 1 and len(config[arg]) > 1:
            config[arg] = [max(config[arg])]
            solver.solve(config, config.keys())
            logger.debug("-" * 72)
            print_config(config)
            domain_size = solver.domain_size(config)

    if not math.isfinite(domain_size):
        raise ConfigurationError("could not determine configuration")

    maximize("inner_message")
    maximize("inner_ecc")
    maximize("outer_ecc")

    if domain_size == 1:
        assert solver.check_constraints(config)
        print_config_args(config)
    else:
        raise ConfigurationError("could not determine configuration")

    return instantiate_config(config)


re_algo = re.compile(r"(RSi16|CIRC)\((\w+=)?[\d+\-*/]+(,\s*(\w+=)?[\d+\-*/]+)*\)")


def parse_algo(algo):
    assert re_algo.match(algo), f"Invalid format: {algo}"
    return eval("ffrs" + algo)


def set_log_verbosity(args):
    verbosity = args.verbosity.get()
    levels = {
        "critical": logging.CRITICAL,
        "warning": logging.WARNING,
        "info": logging.INFO,
        "debug": logging.DEBUG,
    }
    log_level_diff = 0
    par_sub_loggers = {l.name[len(ffrs.par.logger.name) + 1 :]: l for l in ffrs.par.logger.getChildren()}
    for arg in verbosity:
        if arg is None:
            # -v -v
            log_level_diff -= 10
        elif re.match(r"v+$", arg):
            # -vv
            log_level_diff -= 10 * (1 + len(arg))
        elif arg in levels:
            # -vwarning
            log_level_diff = levels[arg] - logger.level
        else:
            # -vsolver=10
            arg_split = arg.split("=")
            if len(arg_split) != 2:
                logger.warning("could not parse verbosity setting '%s'", arg)
                continue

            logger_name, sub_logger_level = arg_split

            if sub_logger_level in levels:
                sub_logger_level = levels[sub_logger_level]
            else:
                try:
                    sub_logger_level = int(sub_logger_level)
                except ValueError:
                    logger.warning("could not parse verbosity setting '%s'", arg)
                    continue

            if logger_name not in par_sub_loggers:
                breakpoint()
                logger.warning("unknown logger '%s'", arg)
                continue

            par_sub_loggers[logger_name].setLevel(sub_logger_level)

    ffrs.par.logger.setLevel(
        max(logging.DEBUG, ffrs.par.logger.level + min(0, logger.level - logging.DEBUG + log_level_diff))
    )
    logger.setLevel(max(logging.DEBUG, logger.level + log_level_diff))


def main():
    cli = CLI()
    args = cli.parse_args()
    set_log_verbosity(args)
    logger.debug("args: %s", args)

    for filename in args.input_file.get():
        if not os.path.isfile(filename):
            logger.critical("input file '%s' does not exist", filename)
            return 1

        output_filename = filename + ".ffrs"

        if os.path.isfile(output_filename):
            if args.force.get():
                logger.info("overwriting '%s'", output_filename)
            else:
                logger.critical("output file '%s' already exists", output_filename)
                return 1

        input_stat = os.stat(filename)

        if input_stat.st_size == 0:
            logger.debug("skipping empty file: '%s'", filename)
            continue

        rs = maximize_interleaving(args, input_stat.st_size)
        logger.info("codec: %s", rs)
        logger.info("buffer size: %s", format_size(rs.message_size))

        # TODO:
        assert rs.message_size >= input_stat.st_size, "Chunking not implemented yet"

        sha1 = hashlib.sha1()
        buffer = ffrs.create_buffer(rs.message_size)
        assert len(buffer) == rs.message_size

        with open(output_filename, "wb") as f:
            f.write(b"\x89ffrs\r\n")
            f.write(rs.serialize())
            f.write(f"\ns:{ffrs.util.b64hex(input_stat.st_size)}\n".encode("ascii"))
            f.write(f"m:{ffrs.util.b64hex(input_stat.st_mtime_ns)}\n".encode("ascii"))
            f.write(b"h:")
            hash_offset = f.tell()
            f.write(base64.urlsafe_b64encode(bytes(sha1.digest_size)).strip(b"="))
            f.write(b"\n\n")

            with open(filename, "rb") as input_file:
                while read := input_file.readinto(buffer):
                    if read < rs.message_size:
                        buffer[read:] = bytes(rs.message_size - read)

                    sha1.update(buffer)
                    encoded = rs.encode(buffer)
                    f.write(encoded)
                f.seek(hash_offset)
                f.write(base64.urlsafe_b64encode(sha1.digest()).strip(b"="))
    logger.info("done")
    return 0


if __name__ == "__main__":
    try:
        quit(main())
    except ffrs.par.FfrsParException as e:
        logger.critical("%s", e)
    except Exception:
        logger.exception("unexpected error")
        quit(1)
