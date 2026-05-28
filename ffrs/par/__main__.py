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

import copy
import itertools
import logging
import math
import re
import sys

import ffrs
import ffrs.par

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
    logger.info("config: {")
    for k in sorted(config):
        val = config[k]
        if val is None:
            val = "\033[35m[[[...]]]\033[0m"
        elif len(val) > 32:
            val = f"\033[35m[[[{len(val)}]]]\033[0m"
        else:
            val = "[" + ", ".join((map(format_size, sorted(val)))) + "]"

        logger.info(f"  {repr(k):23}: {val},")
    logger.info("}\n")


def print_config_args(config):
    assert (
        len(config["inner_block"])
        == len(config["inner_ecc"])
        == len(config["outer_block"])
        == len(config["outer_ecc"])
        == len(config["outer_interleave"])
        == 1
    )
    logger.info(
        f"CIRC16({(*config['inner_block'],)[0]}, {(*config['inner_ecc'],)[0]},"
        f" {(*config['outer_block'],)[0]}, {(*config['outer_ecc'],)[0]}, {(*config['outer_interleave'],)[0]})"
    )


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
        outer_interleave=args.outer_interleave.get(),
        inner_interleaved_ecc=args.outer_interleaved_ecc.get(),
        outer_interleaved_ecc=args.outer_interleaved_ecc.get(),
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
    logger.info("domain size: %s", domain_size)

    if not math.isfinite(domain_size):
        raise ConfigurationError("could not determine configuration")

    if not args.inner_message.is_set() and domain_size > 1 and len(config["inner_message"]) > 1:
        config["inner_message"] = [max(config["inner_message"])]
        solver.solve(config, config.keys())
        logger.info("-" * 72)
        print_config(config)

    if not args.outer_ecc.is_set() and domain_size > 1 and len(config["outer_ecc"]) > 1:
        config["outer_ecc"] = [max(config["outer_ecc"])]
        solver.solve(config, config.keys())
        logger.info("-" * 72)
        print_config(config)


    if domain_size == 1:
        assert solver.check_constraints(config)
        print_config_args(config)
    else:
        logger.critical("* could not determine args")


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

    ffrs.par.logger.setLevel(max(logging.DEBUG, ffrs.par.logger.level + logger.level + log_level_diff))
    logger.setLevel(max(logging.DEBUG, logger.level + log_level_diff))


def main():
    cli = CLI()
    args = cli.parse_args()
    set_log_verbosity(args)
    logger.debug("args: %s", args)

    circ(args)

    return 0


if __name__ == "__main__":
    try:
        quit(main())
    except ffrs.par.FfrsParException as e:
        logger.critical("%s", e)
    except Exception:
        logger.exception("unexpected error")
        quit(1)
