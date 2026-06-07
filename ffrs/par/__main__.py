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

"""File parity tools"""

import math
import re

import ffrs
import ffrs.par
import ffrs.par.exc
import ffrs.util

from . import constraints, cli
from .cli import CLI, log
from .solver import Solver

import ffrs.par.backup
import ffrs.par.benchmark
import ffrs.par.repair


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

    ffrs.par.log.debug(ffrs.util.format_config(config))

    solver = Solver(config.keys(), constraints, free_variables)

    try:
        solver.solve(config)
    except ValueError:
        ffrs.par.log.exception("could not determine configuration")
        return

    ffrs.par.log.debug(ffrs.util.format_config(config))

    domain_size = solver.domain_size(config)
    ffrs.par.log.debug("domain size: %s", domain_size)

    if not math.isfinite(domain_size):
        raise ffrs.par.exc.OptimizationError("could not determine configuration")

    solver.maximize(config, "inner_message")
    solver.maximize(config, "inner_ecc")
    solver.maximize(config, "outer_ecc")

    if domain_size == 1:
        assert solver.check_constraints(config)
        ffrs.util.print_config_args(config)
    else:
        raise ffrs.par.exc.OptimizationError("could not determine configuration")

    return ffrs.util.instantiate_config(config)


re_algo = re.compile(r"(RSi16|CIRC)\((\w+=)?[\d+\-*/]+(,\s*(\w+=)?[\d+\-*/]+)*\)")


def parse_algo(algo):
    assert re_algo.match(algo), f"Invalid format: {algo}"
    return eval("ffrs" + algo)


def main():
    parser = cli.new_parser(__package__, __doc__)

    modules = {
        "backup": ffrs.par.backup,
        "benchmark": ffrs.par.benchmark,
        "repair": ffrs.par.repair,
    }

    subparsers = parser.add_subparsers(dest="command", required=True)
    for name, module in modules.items():
        subparser = subparsers.add_parser(
            name,
            help=module.__doc__.lstrip().split("\n", 1)[0],
            description=module.__doc__,
            formatter_class=parser.formatter_class,
        )
        module.arg_parser(subparser)

    args = parser.parse_args()
    return args.cli_main(args)


if __name__ == "__main__":
    try:
        quit(main())
    except ffrs.par.FfrsParException as e:
        log.critical("%s", e)
    except Exception:
        log.exception("unexpected error")
        quit(1)
