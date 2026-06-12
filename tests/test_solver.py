#  test_solver.py
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
import pytest

import ffrs
import ffrs.par.solver
from ffrs.util import format_config
from ffrs.par.opt import CONSTRAINTS

test_config = dict(
    block=None,
    ecc=None,
    inner_block=range(2, 65536 + 1, 2),
    inner_ecc=set(2**i for i in range(1, 16)),
    inner_interleaved_ecc=None,
    inner_message=range(2, 65536 + 1, 2),
    message=None,
    outer_block=range(2, 65536 + 1, 2),
    outer_ecc_ratio_den=range(1, 1024 + 1),
    outer_ecc_ratio_num=[1],
    outer_ecc=set(2**i for i in range(1, 16)),
    outer_interleave=None,
    outer_interleaved_ecc=None,
    outer_message=range(2, 65536 + 1, 2),
    interleave=None,
)

set_values = dict(inner_message=[4096 - 8, 4096], inner_ecc=[8], outer_ecc=[8, 16, 32, 64, 128, 256])

test_config.update(set_values)

test_config_free = set(test_config.keys()) - set(set_values.keys())


class TestSolver:
    @pytest.mark.parametrize(
        "vars",
        [
            dict(
                message=[2**30],
                outer_ecc_ratio_den=[16],
                interleave=[64],
            ),
            dict(
                block=[2**30],
                outer_ecc_ratio_den=[15],
                interleave=[64],
            ),
        ],
    )
    def test_full_solution(self, vars):
        config = copy.deepcopy(test_config)
        config.update(vars)
        solver = ffrs.par.solver.Solver(config.keys(), CONSTRAINTS, test_config_free | set(vars.keys()))
        solver.solve(config)
        print(format_config(config))
        assert solver.domain_size(config) == 1
        assert solver.check_constraints(config)


@pytest.mark.parametrize(
    "vars",
    [
        dict(
            message=[2**20],
            outer_ecc_ratio_den=[16],
            # interleave=[64],
        ),
        dict(
            block=[2**20],
            outer_ecc_ratio_den=[15],
            # interleave=[64],
        ),
        dict(
            message=[2**30],
            outer_ecc_ratio_den=[16],
            # interleave=[64],
        ),
        dict(
            block=[2**30],
            outer_ecc_ratio_den=[15],
            # interleave=[64],
        ),
    ],
)
class TestSolverPartialConfig:
    def test_partial_solution(self, vars):
        config = copy.deepcopy(test_config)
        config.update(vars)
        solver = ffrs.par.solver.Solver(config.keys(), CONSTRAINTS, test_config_free | set(vars.keys()))
        solver.solve(config)
        print(format_config(config))
        assert 1 < solver.domain_size(config) <= 20736

    def test_partial_solution_continue(self, vars):
        config = copy.deepcopy(test_config)
        config.update(vars)
        solver = ffrs.par.solver.Solver(config.keys(), CONSTRAINTS, test_config_free | set(vars.keys()))
        solver.solve(config)
        print(format_config(config))
        assert 1 < solver.domain_size(config) <= 20736
        assert len(config["outer_ecc"]) > 1

        config["outer_ecc"] = [max(config["outer_ecc"])]
        # solver.free_variables.remove("outer_ecc")
        solver.solve(config, ["outer_ecc"])
        print(format_config(config))
        assert solver.domain_size(config) == 1
        assert solver.check_constraints(config)
