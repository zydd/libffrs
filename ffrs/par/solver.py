#  par/solver.py
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
import math
import re
import time

import ffrs.par
from . import logger as parent_logger

MAX_DOMAIN_SIZE = 65537

logger = parent_logger.getChild("solver")


class SolverError(ffrs.par.FfrsParException):
    pass


class _Equation:
    _re_expression = re.compile(r"[\w\s{}+-]+$")
    _re_var = re.compile(r"\{(\w+)\}")

    @staticmethod
    def solve_for_vars(expr):
        eq = expr.split("==")
        if len(eq) != 2:
            return []

        lhs, rhs = eq
        lhs_match = _Equation._re_expression.match(lhs)
        rhs_match = _Equation._re_expression.match(rhs)
        if not lhs_match and not rhs_match:
            return []

        if lhs_match and rhs_match:
            lhs = _Equation._parse_eq_terms(lhs)
            rhs = _Equation._parse_eq_terms(rhs)

            expr = _Equation._sub(lhs, rhs)
            if not all(re.match(r"\{\w+\}$", term) for sign, term in expr):
                logger.debug("unexpected term: %s", expr)
                return []

            expr = [(sign, term.strip("{}")) for sign, term in expr]
            variables = set(term for sign, term in expr)
            return _Equation._solve_all(variables, expr)

        if not lhs_match:
            lhs, rhs = rhs, lhs

        lhs = _Equation._parse_eq_terms(lhs)

        if not all(re.match(r"\{\w+\}$", term) for sign, term in lhs):
            logger.debug("unexpected term: %s", expr)
            return []

        lhs = [(sign, term.strip("{}")) for sign, term in lhs]

        rhs_vars = set(_Equation._re_var.findall(rhs))

        # Treat rhs as a single term
        rhs = _Equation._re_var.sub(r"\1", rhs)
        rhs = [("+", rhs)]

        expr = _Equation._sub(lhs, rhs)
        variables = set(term for sign, term in lhs)
        return _Equation._solve_all(variables, expr, rhs_vars)

    def _parse_eq_terms(expr):
        expr = re.split(r"([-+])", expr)
        sign = "+"
        res = []
        for term in expr:
            if term in ["+", "-"]:
                sign = term
            else:
                res.append((sign, term.strip()))
        return res

    def _sub(lhs, rhs):
        res = list(lhs)
        for sign, term in rhs:
            if (sign, term) in res:
                res.remove((sign, term))
            else:
                opposite_sign = "+" if sign == "-" else "-"
                res.append((opposite_sign, term))
        return res

    def _solve_for(var, expr):
        res = []
        occurrences = 0
        var_sign = None
        for sign, term in expr:
            if term == var:
                occurrences += 1
                var_sign = sign

        assert occurrences == 1, f"{var} appears more than once in {expr}"

        for sign, term in expr:
            if term == var:
                continue

            if var_sign == "-":
                res.append((sign, term))
            else:
                opposite_sign = "+" if sign == "-" else "-"
                res.append((opposite_sign, term))
        return res

    def _solve_all(vars, expr, other_args=None):
        other_args = other_args or []
        eqs = []
        for var in vars:
            rhs = _Equation._solve_for(var, expr)
            rhs = f"{' '.join(f'{sign} {term}' for sign, term in rhs)}"
            args = [v for v in vars if v != var] + list(other_args)
            func_code = f"lambda {", ".join(args)}: {rhs}"
            eqs.append(
                dict(
                    eq=f"{var} = {rhs}",
                    var=var,
                    func_code=func_code,
                    func=eval(func_code),
                    args=args,
                )
            )
        return eqs


class Solver:
    def __init__(self, variables, constraints, free_variables):
        logger.info("variables: %s", variables)
        logger.info("free variables: %s", free_variables)
        self.free_variables = set(free_variables)
        self._init_constraints(variables, constraints)
        self._init_equivalences()
        self._init_generators()
        self._resolve_input_constraints()

    def _init_constraints(self, variables, constraints):
        self._constraint_evaluations = 0
        self.constraints = constraints
        self._resolved_constraints = set()

        # Constraints for each variable
        self.var_constraints = {}
        for var in variables:
            self.var_constraints[var] = list()
            for i, constr in enumerate(constraints):
                if f"{{{var}}}" in constr:
                    self.var_constraints[var].append(i)

        # Variables needed by each constraint
        self.constraint_vars = [None] * len(constraints)
        for i, constr in enumerate(constraints):
            self.constraint_vars[i] = list()
            for var in variables:
                if f"{{{var}}}" in constr and var not in self.constraint_vars[i]:
                    self.constraint_vars[i].append(var)

        self.constraint_functions = []
        for i, constr in enumerate(constraints):
            args = ", ".join(self.constraint_vars[i])
            constr = constr.format(**dict(zip(self.constraint_vars[i], self.constraint_vars[i])))
            func = f"lambda {args}: {constr}"
            self.constraint_functions.append(eval(func))

    def _init_equivalences(self):
        # Equivalences
        self.equivalences = []
        self.equivalence_functions = []
        self.equivalence_vars = []

        logger.debug("equivalences:")
        for i, constr in enumerate(self.constraints):
            for eq in _Equation.solve_for_vars(constr):
                logger.debug(f"{len(self.equivalences):2} {eq["var"]} = {eq["func_code"]}")
                self.equivalences.append(eq["eq"])
                self.equivalence_functions.append((eq["var"], eq["func"]))
                self.equivalence_vars.append(eq["args"])
        logger.debug("")

    def _init_generators(self):
        # Generators
        re_divisible = re.compile(r" \{(\w+)\} % \{(\w+)\} == 0 $".replace(" ", r"\s*"))
        self.generators = []
        for constr in self.constraints:
            match = re_divisible.match(constr)
            if not match:
                continue

            num = match[1]
            den = match[2]

            self.generators.append((num, den, constr))

    def _resolve_input_constraints(self):
        # Assume that input already satisfies range checks
        re_inequality = re.compile(r"\d+ <= \{\w+\} (<= \d+)$".replace(" ", r"\s*"))
        for i, constr in enumerate(self.constraints):
            if re_inequality.match(constr):
                logger.debug("assume resolved: %2d '%s'", i, constr)
                self._resolved_constraints.add(i)

    def _propagate_equivalences(self, config):
        updated = False
        for i, (var, func) in enumerate(self.equivalence_functions):
            arg_domain_size = math.prod(
                len(config[arg]) if config[arg] is not None else float("inf") for arg in self.equivalence_vars[i]
            )
            var_domain_size = len(config[var]) if config[var] is not None else float("inf")

            logger.info("equivalence %2s: '%s'", i, self.equivalences[i])
            logger.debug("domain size: var: %s arg: %s", var_domain_size, arg_domain_size)

            if not (
                arg_domain_size < var_domain_size
                and arg_domain_size <= MAX_DOMAIN_SIZE
                # and var_domain_size > MAX_DOMAIN_SIZE
                and var in self.free_variables
            ):
                continue

            # self.free_variables.remove(var)
            valid_set = set(
                func(*values) for values in itertools.product(*(config[arg] for arg in self.equivalence_vars[i]))
            )
            logger.info(
                "%-23s %5s -> %-5d     '%s'",
                var,
                var_domain_size,
                len(valid_set),
                self.equivalences[i],
            )
            if config[var] and len(config[var]) <= 32:
                logger.debug("removing: %s", config[var] - valid_set)
                logger.debug("remaining: %s", valid_set)

            config[var] = valid_set

            # Need to re-evaluate constraints for new set
            self._resolve_constraints(config, self.var_constraints[var])
            # self._resolved_constraints -= set(self.var_constraints[var])

            updated = True
        return updated

    def _gen_values(self, config):
        for i, (num, den, gen) in enumerate(self.generators):
            if den not in self.free_variables:
                continue
            if config[num] is None:
                continue
            if config[den] is None:
                continue
            if not all(v & (v - 1) == 0 for v in config[num]):
                continue

            den_min_log = math.log2(min(config[den]))

            den_v_max_log = round(min(math.log2(max(config[den])), max(math.log2(v) for v in config[num])))
            if den_v_max_log > len(config[den]):
                continue

            den_v = {2**i for i in range(math.ceil(den_min_log), round(den_v_max_log) + 1)}

            if len(den_v) < len(config[den]):
                logger.info(f"{den:23} {len(config[den]):5} -> {len(den_v):<5}     '{gen}'")
                config[den] = den_v
                self._resolve_constraints(config, self.var_constraints[den])

    def _eval_constraint(self, config, idx):
        logger.info("constraint: '%s'", self.constraints[idx])
        evaluations = 0
        updated = False
        valid = [set() for _ in range(len(self.constraint_vars[idx]))]
        constraint = self.constraint_functions[idx]
        for var_i, var in enumerate(self.constraint_vars[idx]):
            # logger.warn("var: %s", var)
            constr_domain = [config[v] for v in self.constraint_vars[idx]]
            for val in config[var]:
                # logger.warn("val: %s", var)
                constr_domain[var_i] = [val]
                for constr_args in itertools.product(*constr_domain):
                    try:
                        evaluations += 1
                        constr_eval = constraint(*constr_args)
                    except Exception:
                        logger.exception(
                            "error evaluating constraint '%s' with values %s",
                            self.constraints[idx],
                            dict(zip(self.constraint_vars[idx], constr_args)),
                        )
                        breakpoint()
                        continue
                    if constr_eval is True:
                        valid[var_i].add(val)
                        break

        for i, var in enumerate(self.constraint_vars[idx]):
            # log = logger.debug if len(valid[i]) == len(config[var]) else logger.info
            if len(valid[i]) == len(config[var]):
                continue
            logger.info(f"{var:23} {len(config[var]):5} -> {len(valid[i]):<5}     '{self.constraints[idx]}'")

            if config[var] and len(config[var]) <= 32:
                logger.debug("removing: %s", config[var] - valid[i])
                logger.debug("remaining: %s", valid[i])

            if len(valid[i]) == 0:
                raise ValueError(
                    f"Constraint '{self.constraints[idx]}' is unsatisfiable with current configuration for variable '{var}'"
                )

            config[var] = valid[i]
            updated = True

        self._constraint_evaluations += evaluations
        if evaluations > 10000:
            logger.warning("evaluations: %s constraint: '%s'", evaluations, self.constraints[idx])
        else:
            logger.info("evaluations: %s", evaluations)
        return updated

    def _resolve_constraints(self, config, constraints=None):
        constraints = constraints or (set(range(len(self.constraints))) - self._resolved_constraints)
        # constraints = range(len(self.constraints))  # force re-evaluation

        domain_sizes = map(lambda idx: (self._constraint_domain_size(config, idx), idx), constraints)
        domain_sizes = sorted(filter(lambda x: x[0] <= MAX_DOMAIN_SIZE, domain_sizes))
        updated = 0
        for domain_size, idx in domain_sizes:
            updated += self._eval_constraint(config, idx)
            self._resolved_constraints.add(idx)
        return updated

    def _constraint_domain_size(self, config, idx):
        domain_size = math.prod(
            len(config[var]) if config[var] is not None else float("inf") for var in self.constraint_vars[idx]
        )
        return domain_size

    def observable_domain_size(self, config):
        return math.prod(len(config[var]) if config[var] is not None else 1 for var in config)

    def domain_size(self, config):
        return math.prod(len(config[var]) if config[var] is not None else float("inf") for var in config)

    def solve(self, config, vars=None):
        if vars:
            for var in vars:
                self._resolved_constraints -= set(self.var_constraints[var])

        for k, v in config.items():
            if v and type(v) is not set:
                config[k] = set(v)

        domain_size = self.observable_domain_size(config)
        logger.info("observable domain size: %s", domain_size)

        while domain_size > 1:
            initial_domain_size = domain_size

            while self._propagate_equivalences(config):
                pass
            # self._gen_values(config)
            self._resolve_constraints(config)

            domain_size = self.observable_domain_size(config)
            logger.info("observable domain size: %s", domain_size)
            if domain_size == initial_domain_size:
                break

        logger.info("constraint evaluations: %s", self._constraint_evaluations)
        return domain_size

    def check_constraints(self, config):
        for idx, constr in enumerate(self.constraints):
            domain_size = self._constraint_domain_size(config, idx)
            if domain_size > MAX_DOMAIN_SIZE:
                logger.info(f"constraint '{constr}' is not satisfied with current configuration")
                return False

            for constr_args in itertools.product(*(config[var] for var in self.constraint_vars[idx])):
                constr_eval = self.constraint_functions[idx](*constr_args)
                if not constr_eval:
                    logger.info(f"constraint '{constr}' is not satisfied with current configuration")
                    return False
        else:
            return True

    def branch_config(self, config, arg, value):
        t0 = time.time()
        # logger.debug("branch %s = %s", arg, value)
        # config2 = copy.deepcopy(config)  # slower than making a new dict
        config2 = {k: set(v) if v is not None else None for k, v in config.items()}
        config2[arg] = [value]
        self.solve(config2, [arg])
        logger.debug("branching %s = %s took %.2f ms", arg, value, (time.time() - t0) * 1000)
        return config2

    def _optimize(self, mode, config, arg):
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
                    config2 = self.branch_config(config, arg, val)
                    logger.debug("optimizing %s = %s (index: %d) took %.2f ms", arg, val, i, (time.time() - t0) * 1000)
                    return config2
                except ValueError:
                    pass
            raise SolverError(f"could not optimize {arg}")
        return config

    def minimize(self, config, arg):
        logger.info("minimize %s", arg)
        return self._optimize("min", config, arg)

    def maximize(self, config, arg):
        logger.info("maximize %s", arg)
        return self._optimize("max", config, arg)
