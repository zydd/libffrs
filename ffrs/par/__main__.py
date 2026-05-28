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
import re

import ffrs

from . import constraints, print_config, print_config_args, logger
from .cli import CLI
from .solver import Solver


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
        interleaved_ecc=args.interleaved_ecc.get(),
    )

    free_variables = {k for k in config if k not in args or not getattr(args, k).has_value()}
    if not not args.ecc_ratio.has_value():
        free_variables.update("outer_ecc_ratio_num", "outer_ecc_ratio_den")

    print_config(config)

    solver = Solver(config.keys(), constraints, free_variables)

    solver.solve(config)
    print(solver._iterations)
    print_config(config)

    if not args.outer_ecc.is_set() and solver.domain_size(config) > 1 and len(config["outer_ecc"]) > 1:
        for ecc in config["outer_ecc"]:
            config_ecc = copy.deepcopy(config)
            config_ecc["outer_ecc"] = [ecc]
            solver.solve(config_ecc, config_ecc.keys())
            print_config_args(config_ecc)

    # if solver.check_constraints(config) and 0 < solver.domain_size(config) < 10:
    #     configs = list(
    #         itertools.product(
    #             config["inner_block"],
    #             config["inner_ecc"],
    #             config["outer_block"],
    #             config["outer_ecc"],
    #             config["outer_interleave"],
    #         )
    #     )
    #     for cfg in configs:
    #         print(cfg)


re_algo = re.compile(r"(RSi16|CIRC)\((\w+=)?[\d+\-*/]+(,\s*(\w+=)?[\d+\-*/]+)*\)")


def parse_algo(algo):
    assert re_algo.match(algo), f"Invalid format: {algo}"
    return eval("ffrs" + algo)


def main():
    cli = CLI()
    args = cli.parse_args()
    print(args)

    verbosity = args.verbosity.get()

    levels = {
        "critical": logging.CRITICAL,
        "warning": logging.WARNING,
        "info": logging.INFO,
        "debug": logging.DEBUG,
    }
    level = logging.WARNING
    sub_loggers = {l.name[len(logger.name) + 1 :]: l for l in logger.getChildren()}
    for arg in verbosity:
        if arg is None:
            # -v -v
            level = -10
        elif re.match(r"v+$", arg):
            # -vv
            level = -10 * len(arg)
        elif arg in levels:
            # -vwarning
            level = levels[arg]
        else:
            # -vsolver=10
            arg_split = arg.split("=")
            if len(arg_split) != 2:
                logger.warning("Could not parse verbosity setting '%s'", arg)
                continue

            logger_name, sub_logger_level = arg_split

            if sub_logger_level in levels:
                sub_logger_level = levels[sub_logger_level]
            else:
                try:
                    sub_logger_level = int(sub_logger_level)
                except ValueError:
                    logger.warning("Could not parse verbosity setting '%s'", arg)
                    continue

            if logger_name not in sub_loggers:
                logger.warning("Unknown logger '%s'", arg)
                continue

            sub_loggers[logger_name].setLevel(sub_logger_level)

    print("verbosity:", type(verbosity), level)
    logger.setLevel(max(logging.DEBUG, level))

    circ(args)


if __name__ == "__main__":
    main()
