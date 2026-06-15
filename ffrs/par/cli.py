#
#  par/cli.py
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

import argparse
import inspect
import logging
import re
import sys

import ffrs.par


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


log = logging.getLogger("__main__")
log.setLevel(logging.INFO)

stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setLevel(logging.DEBUG)
stdout_handler.addFilter(MaxLevelFilter(logging.WARNING))
stdout_handler.setFormatter(ColorFormatter())
log.addHandler(stdout_handler)

stderr_handler = logging.StreamHandler(sys.stderr)
stderr_handler.setLevel(logging.WARNING)
stderr_handler.setFormatter(ColorFormatter())
log.addHandler(stderr_handler)


_UNSET = object()
_DEFAULT_LIB_LOG_LEVEL = ffrs.par.log.level
_DEFAULT_CLI_LOG_LEVEL = log.level


class Value:
    _is_default = False

    def __init__(self, value, value_parser=None):
        self.value_text = value
        self._value = value_parser(value) if value_parser else value

    def __repr__(self):
        return f"({"default" if self._is_default else "value"}: {repr(self._value)})"

    def get(self, default=_UNSET):
        if self._value is not None:
            return self._value
        else:
            if default is _UNSET:
                raise ValueError("no value set")
            return default

    def is_default(self):
        return self._is_default

    def is_set(self):
        return self.has_value() and not self.is_default()

    def has_value(self):
        return self._value is not None

    def append(self, val):
        if self._is_default:
            self._is_default = False
            self._value = [val]
        else:
            self._value.append(val)


class DEFAULT(Value):
    _is_default = True


class HelpFormatter(argparse.HelpFormatter):
    def __init__(self, prog, indent_increment=2, max_help_position=35, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)

    def _get_help_string(self, action):
        help = action.help
        if action.default is not argparse.SUPPRESS and action.default is not None:
            default = action.default
            if isinstance(default, DEFAULT):
                default = default.value_text
            help += f" [default: {default}]"
        if action.choices and action.metavar:
            help += f" [choices: {', '.join(str(c) for c in action.choices)})]"
        elif action.type and hasattr(action.type, "choices"):
            help += f" [choices: {', '.join(str(c) for c in sorted(action.type.choices))}]"

        return help


class UpdateAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        current = getattr(namespace, self.dest, None)
        if isinstance(current, DEFAULT):
            current = {}
        current.update(values)
        setattr(namespace, self.dest, current)


class CLI:
    DEFAULT = globals()["DEFAULT"]

    def __init__(self, parser: argparse.ArgumentParser):
        self.checks = []
        self.parser = parser
        self.cli_verbosity = _DEFAULT_CLI_LOG_LEVEL
        self.lib_verbosity = _DEFAULT_LIB_LOG_LEVEL

    def check(self, func, *args):
        self.checks.append((func, args))

    def warn(self, namespace, func, message):
        if not func(namespace):
            print(message, file=sys.stderr)

    def check_either(self, namespace, *args):
        dests = [arg.lstrip("-").replace("-", "_") for arg in args]
        if not any(not getattr(namespace, arg).is_default() for arg in dests):
            raise ValueError(f"At least one of {', '.join(args)} must be specified")

    def check_mutually_exclusive(self, namespace, *args):
        dests = [arg.lstrip("-").replace("-", "_") for arg in args]
        count = sum(getattr(namespace, arg).is_set() for arg in dests)
        if count > 1:
            raise ValueError(f"Only one of {', '.join(args)} can be specified")

        if count == 1:
            for arg in dests:
                if getattr(namespace, arg).is_default():
                    setattr(namespace, arg, Value(None))

    def run_checks(self, namespace):
        for func, args in self.checks:
            func(namespace, *args)

    def parse_args(self, argv=None):
        args = self.parser.parse_args(argv)

        for key, value in args.__dict__.items():
            if not isinstance(value, Value):
                setattr(args, key, Value(value))

        try:
            self.run_checks(args)
        except ValueError as e:
            # import traceback; traceback.print_exc()
            self.parser.error(str(e))

        return args

    LOG_LEVELS = {
        "critical": logging.CRITICAL,
        "warning": logging.WARNING,
        "info": logging.INFO,
        "debug": logging.DEBUG,
    }

    def parse_verbosity(self, arg):
        log_level_diff = 0
        par_sub_logs = {l.name[len(ffrs.par.log.name) + 1 :]: l for l in ffrs.par.log.getChildren()}

        if arg == "":
            # -v -v
            log_level_diff -= 10
        elif re.match(r"v+$", arg):
            # -vv
            log_level_diff -= 10 * (1 + len(arg))
        elif arg in CLI.LOG_LEVELS:
            # -vwarning
            log_level_diff = CLI.LOG_LEVELS[arg] - _DEFAULT_CLI_LOG_LEVEL
        else:
            # -vsolver=10
            arg_split = arg.split("=")
            if len(arg_split) != 2:
                log.warning("could not parse verbosity setting '%s'", arg)
                return

            log_name, sub_log_level = arg_split

            if sub_log_level in CLI.LOG_LEVELS:
                sub_log_level = CLI.LOG_LEVELS[sub_log_level]
            else:
                try:
                    sub_log_level = int(sub_log_level)
                except ValueError:
                    log.warning("could not parse verbosity setting '%s'", arg)
                    return

            if log_name not in par_sub_logs:
                log.warning("unknown log '%s'", arg)
                return

            par_sub_logs[log_name].setLevel(sub_log_level)

        self.lib_verbosity = max(
            logging.DEBUG, self.lib_verbosity + min(0, self.cli_verbosity - logging.DEBUG + log_level_diff)
        )
        self.cli_verbosity = max(logging.DEBUG, self.cli_verbosity + log_level_diff)

        ffrs.par.log.setLevel(self.lib_verbosity)
        log.setLevel(max(logging.DEBUG, self.cli_verbosity))
        ffrs.set_logger(ffrs.par.log.getChild("rs"))
        return arg


def parse_size(size_str: str) -> int:
    size_str = size_str.strip().lower()

    multipliers = {
        "k": 1024,
        "m": 1024**2,
        "g": 1024**3,
        "t": 1024**4,
    }

    total = 0
    for chunk in re.split(r"(?<=[a-zA-Z])(?=\d)", size_str):
        if chunk[-1] in multipliers:
            number = float(chunk[:-1]) if "" in chunk[:-1] else int(chunk[:-1])
            total += int(number * multipliers[chunk[-1]])
        else:
            total += int(chunk)
    return total


def parse_ratio(ratio_str: str) -> tuple[float, float]:
    ratio_str = ratio_str.strip().lower()

    if "/" in ratio_str:
        parts = ratio_str.split("/")
    elif ":" in ratio_str:
        parts = ratio_str.split(":")
    else:
        raise ValueError(f"Invalid ratio format: {ratio_str}")

    if len(parts) != 2:
        raise ValueError(f"Invalid ratio format: {ratio_str}")

    numerator = parse_size(parts[0])
    denominator = parse_size(parts[1])

    return numerator, denominator


def parse_choice_list(choices):
    choices = set(choices)

    def parser(value):
        items = value.split(",")

        invalid = [x for x in items if x not in choices]
        if invalid:
            raise argparse.ArgumentTypeError(
                f"invalid choice(s): {', '.join(invalid)} " f"(choose from {', '.join(sorted(choices))})"
            )

        return items

    parser.choices = choices
    return parser


def parse_int_list(value):
    return [int(i) for i in value.split(",")]


def parse_size_list(value):
    if value == "*":
        return None
    return [parse_size(i) for i in value.split(",")]


def parse_ratio_list(value):
    return [parse_ratio(i) for i in value.split(",")]


_re_algo = re.compile(r"(RSi16|CIRC)\((\w+=)?[\d+\-*/]+(,\s*(\w+=)?[\d+\-*/]+)*\)")


def parse_algo(algo):
    if not _re_algo.match(algo):
        raise ValueError(f"Invalid format: {algo}")
    return eval("ffrs" + algo)


class OuterScopeGetter(object):
    NIL = object()

    def __init__(self):
        frame = inspect.currentframe()
        if frame is None or frame.f_back is None or frame.f_back.f_back is None:
            raise RuntimeError("cannot inspect stack frames")
        self.frame = frame.f_back.f_back

    def __getattr__(self, name):
        # frame = self.frame
        # while frame is not None:

        value = self.frame.f_locals.get(name, self.NIL)
        if value is self.NIL:
            value = self.frame.f_globals.get(name, self.NIL)
        if value is not self.NIL:
            return value

        # frame = frame.f_back
        # raise AttributeError(repr(name) + " not found in outer scope")


def create_parser(prog, description, parser=None):
    parser = parser or argparse.ArgumentParser(prog=prog, description=description, formatter_class=HelpFormatter)
    return CLI(parser)


def cli_main(parser, main_function, argv):
    try:
        ffrs.set_logger(ffrs.par.log.getChild("rs"))

        args = parser.parse_args(argv)

        log.debug(args)

        return main_function(args)

    except ffrs.par.FfrsParException as e:
        log.critical("%s", e)
        return 1
    except Exception:
        log.exception("unexpected error")
        return 2


def main():
    outer = OuterScopeGetter()
    parser = outer.create_parser()
    quit(cli_main(parser, outer.main, sys.argv[1:]))
