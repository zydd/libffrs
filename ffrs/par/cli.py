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
_DEFAULT_MOD_LOG_LEVEL = ffrs.par.log.level
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


class CLI:
    DEFAULT = globals()["DEFAULT"]

    def __init__(self, parser, main_function):
        self.checks = []
        self.parser = parser
        self.main_function = main_function

        self.parser.set_defaults(cli_main=self.main)

        self.parser.add_argument(
            "-v",
            "--verbosity",
            metavar="level",
            action="append",
            nargs="?",
            default=DEFAULT("info", lambda a: [a]),
            help="Increase/set verbosity level",
        )

        codec_args = self.parser.add_argument_group("codec")
        codec_args.add_argument(
            "--codec",
            action="store",
            default="circ16",
            choices=["rsi16", "circ16"],
            help="Method use to compute parity",
        )

        self.check(lambda args: args.codec != "rsi16" or self.check_either(args, "--block", "--message"))
        self.check(self.check_either, "--block", "--message")
        # self.check(self.check_mutually_exclusive, "--block", "--message")
        codec_args.add_argument(
            "-b",
            "--block",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Block size to use when encoding. block = message + ecc",
        )
        codec_args.add_argument(
            "-m",
            "--message",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Payload size in a parity block",
        )

        self.check(self.check_mutually_exclusive, "--ecc", "--ecc-ratio")
        codec_args.add_argument(
            "-e",
            "--ecc",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="ECC size in a parity block",
        )
        codec_args.add_argument(
            "-r",
            "--ecc-ratio",
            metavar="ratio",
            action="store",
            default=DEFAULT("1/15,1/16", parse_ratio_list),
            type=parse_ratio_list,
            help="ECC size as a fraction of the message size",
        )

        codec_args.add_argument(
            "--interleave",
            metavar="size",
            action="store",
            # default=DEFAULT("1", parse_size_list),
            type=parse_size_list,
            help="Interleave multiple blocks",
        )

        codec_args.add_argument(
            "--primitive",
            metavar="int",
            action="store",
            default=DEFAULT(3),
            type=int,
            help="Primitive n-th root of unity used to generate GF(65537)",
        )

        self.arg_circ()
        self.arg_simd()

    def arg_simd(self):
        simd = self.parser.add_argument_group("simd")
        simd_choices = ["x4", "x8", "x16", "sse", "avx2", "avx512"]
        self.check(self.check_mutually_exclusive, "--simd", "--no-simd")
        simd.add_argument(
            "--simd",
            metavar="simd",
            nargs="?",
            const="auto",
            default=DEFAULT("auto"),
            type=parse_choice_list(["auto"] + simd_choices),
            help="Enable SIMD optimizations",
        )
        simd.add_argument(
            "--no-simd",
            metavar="simd",
            nargs="?",
            const="*",
            type=parse_choice_list(["*"] + simd_choices),
            help="Disable SIMD optimizations",
        )

    def arg_circ(self):
        circ = self.parser.add_argument_group("circ16")

        def add_circ_arg(*args, **kwargs):
            if "dest" in kwargs:
                dest = kwargs["dest"]
            else:
                assert args[-1].startswith("--")
                dest = args[-1][2:].replace("-", "_")

            self.check(
                self.warn,
                lambda args: args.codec.get() == "circ16" or getattr(args, dest).is_default(),
                "circ16 options have no effect when codec != circ16",
            )
            circ.add_argument(*args, **kwargs)

        add_circ_arg(
            "--inner-block",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Block size used by the inner codec",
        )
        add_circ_arg(
            "--inner-message",
            metavar="size",
            action="store",
            default=DEFAULT("508,512,1020,1024,2040,2048,4088,4096", parse_size_list),
            type=parse_size_list,
            help="Payload size used by the inner codec",
        )
        add_circ_arg(
            "--inner-ecc",
            metavar="size",
            action="store",
            default=DEFAULT("4,8", parse_size_list),
            type=parse_size_list,
            help="ECC size of the inner codec",
        )
        add_circ_arg(
            "--inner-interleaved-ecc",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Total ECC size of the inner codec (inner_ecc * outer_block * interleave)",
        )

        add_circ_arg(
            "--outer-block",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Block size used by the outer codec",
        )
        add_circ_arg(
            "--outer-message",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Payload size used by the outer codec",
        )
        add_circ_arg(
            "--outer-ecc",
            metavar="size",
            action="store",
            default=DEFAULT("16,32,64,128,256", parse_size_list),
            type=parse_size_list,
            help="ECC size of the outer codec",
        )
        add_circ_arg(
            "--outer-interleaved-ecc",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Total ECC size of the outer codec (outer_ecc * inner_message * interleave)",
        )
        add_circ_arg(
            "--outer-interleave",
            metavar="count",
            action="store",
            # default=DEFAULT("1", parse_int_list),
            type=parse_int_list,
            help="Number of interleaved blocks of the outer codec",
        )

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

    def main(self, args=None):
        if args is None:
            args = self.parser.parse_args()

        for key, value in args.__dict__.items():
            if not isinstance(value, Value):
                setattr(args, key, Value(value))

        try:
            self.run_checks(args)
        except ValueError as e:
            # import traceback; traceback.print_exc()
            self.parser.error(str(e))

        set_log_verbosity(args)
        log.debug(args)

        return self.main_function(args)


def new_parser(prog, desc):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=desc,
        formatter_class=HelpFormatter,
    )
    return parser


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


def set_log_verbosity(args):
    verbosity = args.verbosity.get()
    levels = {
        "critical": logging.CRITICAL,
        "warning": logging.WARNING,
        "info": logging.INFO,
        "debug": logging.DEBUG,
    }
    log_level_diff = 0
    par_sub_logs = {l.name[len(ffrs.par.log.name) + 1 :]: l for l in ffrs.par.log.getChildren()}
    for arg in verbosity:
        if arg is None:
            # -v -v
            log_level_diff -= 10
        elif re.match(r"v+$", arg):
            # -vv
            log_level_diff -= 10 * (1 + len(arg))
        elif arg in levels:
            # -vwarning
            log_level_diff = levels[arg] - _DEFAULT_CLI_LOG_LEVEL
        else:
            # -vsolver=10
            arg_split = arg.split("=")
            if len(arg_split) != 2:
                log.warning("could not parse verbosity setting '%s'", arg)
                continue

            log_name, sub_log_level = arg_split

            if sub_log_level in levels:
                sub_log_level = levels[sub_log_level]
            else:
                try:
                    sub_log_level = int(sub_log_level)
                except ValueError:
                    log.warning("could not parse verbosity setting '%s'", arg)
                    continue

            if log_name not in par_sub_logs:
                log.warning("unknown log '%s'", arg)
                continue

            par_sub_logs[log_name].setLevel(sub_log_level)

    ffrs.par.log.setLevel(
        max(logging.DEBUG, _DEFAULT_MOD_LOG_LEVEL + min(0, _DEFAULT_CLI_LOG_LEVEL - logging.DEBUG + log_level_diff))
    )
    log.setLevel(max(logging.DEBUG, _DEFAULT_CLI_LOG_LEVEL + log_level_diff))


class OuterScopeGetter(object):
    def __getattribute__(self, name):
        frame = inspect.currentframe()
        if frame is None:
            raise RuntimeError("cannot inspect stack frames")
        sentinel = object()
        frame = frame.f_back.f_back
        if frame is not None:
            value = frame.f_locals.get(name, sentinel)
            if value is sentinel:
                value = frame.f_globals.get(name, sentinel)
            if value is not sentinel:
                return value
            # frame = frame.f_back
        raise AttributeError(repr(name) + " not found in outer scope")


def main():
    outer = OuterScopeGetter()
    try:
        parser = new_parser(outer.__loader__.name, outer.__doc__)
        cli_parser = outer.arg_parser(parser)
        quit(cli_parser.main())
    except ffrs.par.FfrsParException as e:
        log.critical("%s", e)
    except Exception:
        log.exception("unexpected error")
        quit(1)
