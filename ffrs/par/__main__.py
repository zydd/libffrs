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

import itertools
import re

import ffrs

from . import constraints, print_config
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

    # if not args.outer_ecc.is_set() and solver.domain_size(config) > 1 and len(config["outer_ecc"]) > 1:
    #     config["outer_ecc"] = [max(config["outer_ecc"])]
    #     solver.solve(config, config.keys())
    #     print_config()

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
    circ(args)


if __name__ == "__main__":
    main()
