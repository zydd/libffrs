#  par/opt.py
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

"""Full backup + optional extra parity"""

import math
import time

import ffrs
import ffrs.par.exc
import ffrs.util

from .solver import Solver
from . import log as parent_log

log = parent_log.getChild("opt")


CONSTRAINTS = [
    "2 <= {block}",
    "2 <= {message}",
    "2 <= {ecc}",
    #
    "2 <= {inner_block} <= 65536",
    "2 <= {inner_message} <= 4096",
    "2 <= {inner_ecc} <= 32768",
    #
    "2 <= {outer_block} <= 65536",
    "2 <= {outer_message} <= 65536",
    "2 <= {outer_ecc} <= 32768",
    "1 <= {interleave}",
    "1 <= {outer_interleave}",
    #
    "127 <= {inner_message} // {inner_ecc} <= 512",
    #
    "{block} == {message} + {ecc}",
    "{inner_block} == {inner_message} + {inner_ecc}",
    "{outer_block} == {outer_message} + {outer_ecc}",
    #
    "{block} == {inner_block} * {outer_block} * {interleave}",
    "{message} == {outer_message} * {outer_interleave}",
    "{outer_interleaved_ecc} == {outer_ecc} * {outer_interleave}",
    "{inner_interleaved_ecc} == {outer_block} * {inner_ecc} * {interleave}",
    "{ecc} == {outer_interleaved_ecc} + {inner_interleaved_ecc}",
    #
    "{outer_message} * {outer_ecc_ratio_num} == {outer_ecc} * {outer_ecc_ratio_den}",
    "{message} * {outer_ecc_ratio_num} == {outer_interleaved_ecc} * {outer_ecc_ratio_den}",
    "{outer_block} * {outer_ecc_ratio_num} == {outer_ecc} * ({outer_ecc_ratio_den} + 1)",
    #
    "{outer_ecc} & ({outer_ecc} - 1) == 0",
    #
    "{message} % {inner_message} == 0",
    "{message} % {outer_message} == 0",
    "{message} % {interleave} == 0",
    "{message} % {outer_interleave} == 0",
    "{message} % ({inner_message} * {outer_message}) == 0",
    # "{message} % {ecc} == 0",  # not applicable for CIRC
    "{message} % {outer_interleaved_ecc} == 0",
    "{inner_message} % {inner_ecc} == 0",
    "{outer_message} % {outer_ecc} == 0",
    #
    "{block} % {inner_block} == 0",
    "{block} % {outer_block} == 0",
    "{block} % {interleave} == 0",
    "{block} % ({inner_block} * {outer_block}) == 0",
    #
    "{interleave} == {message} // ({inner_message} * {outer_message})",
    "{interleave} == {block} // ({inner_block} * {outer_block})",
    "{outer_interleave} == {interleave} * {inner_message}",
]


def maximize_interleaving(args, file_size):
    t0 = time.time()
    file_size = max(file_size, 8192)

    config = dict(
        block=args.block.get(None),
        message=None,  # choose based on file_size
        ecc=args.ecc.get(None),
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
        interleave=args.interleave.get(None),
        inner_interleaved_ecc=args.outer_interleaved_ecc.get(None),
        outer_interleaved_ecc=args.outer_interleaved_ecc.get(None),
        outer_interleave=args.outer_interleave.get(None),
    )

    free_variables = {k for k in config if k not in args or not getattr(args, k).has_value()}
    free_variables |= {"message", "inner_message"}
    if args.ecc_ratio.has_value():
        free_variables.remove("outer_ecc_ratio_den")
        free_variables.remove("outer_ecc_ratio_num")

    log.debug(ffrs.util.format_config(config))

    extra_constraints = [
        f"{file_size} <= {{message}} <= {2 * file_size}",
    ]

    solver = Solver(config.keys(), CONSTRAINTS + extra_constraints, free_variables)

    try:
        solver.solve(config)
    except ValueError:
        log.exception("could not determine configuration")
        # return
        raise

    for outer_ecc in sorted(config["outer_ecc"], reverse=True):
        try:
            log.debug("branch outer_ecc")
            config2 = solver.branch_config(config, "outer_ecc", outer_ecc)
            log.debug(ffrs.util.format_config(config2))

            log.debug("compute outer_interleave")
            max_padding = 0.05 if file_size <= 16000 else 0.01
            config2["outer_interleave"] = [
                math.ceil(file_size / outer_message) + i
                for outer_message in config2["outer_message"]
                for i in range(min(100, math.ceil(file_size * max_padding / outer_message)))
            ]
            solver.solve(config2, ["outer_interleave"])
            log.debug(ffrs.util.format_config(config2))

            config2 = solver.maximize(config2, "inner_message")
            log.debug(ffrs.util.format_config(config2))

            config2 = solver.minimize(config2, "inner_ecc")
            log.debug(ffrs.util.format_config(config2))

            config2 = solver.minimize(config2, "message")
            log.debug(ffrs.util.format_config(config2))

            log.info("solver iterations: %s", solver._constraint_evaluations)
            log.info("max message ratio: %s", max(config2["message"]) / file_size)
            log.info("max hash ratio: %s", max(config2["inner_interleaved_ecc"]) / file_size)
            config = config2
            break
        except ValueError:
            log.debug("could not set outer_ecc to %s", outer_ecc)
    else:
        raise ffrs.par.exc.OptimizationError("could not optimize interleaving")

    log.debug(ffrs.util.format_config(config))
    t1 = time.time() - t0
    if t1 > 10e-3:
        log.warning("time to solve: %.2f [ms]", t1 * 1000)
    return ffrs.util.instantiate_config(config)


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

    solver = Solver(config.keys(), CONSTRAINTS, free_variables)

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
