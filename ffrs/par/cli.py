#  par/__init__.py
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
import sys


def parse_size(size_str: str) -> int:
    size_str = size_str.strip().lower()

    multipliers = {
        "k": 1024,
        "m": 1024**2,
        "g": 1024**3,
        "t": 1024**4,
    }

    if size_str[-1] in multipliers:
        number = float(size_str[:-1]) if "" in size_str[:-1] else int(size_str[:-1])
        return int(number * multipliers[size_str[-1]])
    else:
        return int(size_str)


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


class CLI:
    class Value:
        _is_default = False

        def __init__(self, value, value_parser=None):
            self.value_text = value
            self._value = value_parser(value) if value_parser else value

        def __repr__(self):
            return f"({"default" if self._is_default else "value"}: {repr(self._value)})"

        def get(self, default=None):
            return self._value if self._value is not None else default

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
                if isinstance(default, CLI.DEFAULT):
                    default = default.value_text
                help += f" [default: {default}]"
            if action.choices and action.metavar:
                help += f" [choices: {', '.join(str(c) for c in action.choices)})]"
            elif action.type and hasattr(action.type, "choices"):
                help += f" [choices: {', '.join(str(c) for c in sorted(action.type.choices))}]"

            return help

    def __init__(self):
        self.checks = []
        self.main_parser = argparse.ArgumentParser(
            prog="ffrs.par",
            description="File parity",
            formatter_class=self.HelpFormatter,
        )

        self.common_parser = argparse.ArgumentParser(add_help=False)

        self.common_parser.add_argument(
            "-v",
            "--verbosity",
            metavar="level",
            action="append",
            nargs="?",
            default=CLI.DEFAULT("info", lambda a: [a]),
            help="Increase/set verbosity level",
        )

        self.common_parser.add_argument(
            "--codec",
            action="store",
            default="circ16",
            choices=["rsi16", "circ16"],
            help="Method use to compute parity",
        )

        self.check(lambda args: args.codec != "rsi16" or self.check_either(args, "--block", "--message"))
        self.check(self.check_either, "--block", "--message")
        # self.check(self.check_mutually_exclusive, "--block", "--message")
        self.common_parser.add_argument(
            "-b",
            "--block",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Block size to use when encoding. block = message + ecc",
        )
        self.common_parser.add_argument(
            "-m",
            "--message",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Payload size in a parity block",
        )

        self.check(self.check_mutually_exclusive, "--ecc", "--ecc-ratio")
        self.common_parser.add_argument(
            "-e",
            "--ecc",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="ECC size in a parity block",
        )
        self.common_parser.add_argument(
            "-r",
            "--ecc-ratio",
            metavar="ratio",
            action="store",
            # default=CLI.DEFAULT("1/16", parse_ratio_list),
            type=parse_ratio_list,
            help="ECC size as a fraction of the message size",
        )

        self.common_parser.add_argument(
            "--interleave",
            metavar="size",
            action="store",
            # default=CLI.DEFAULT("1", parse_size_list),
            type=parse_size_list,
            help="Interleave multiple blocks",
        )

        self.arg_simd()

        self.common_parser.add_argument(
            "--primitive",
            metavar="int",
            action="store",
            default=CLI.DEFAULT(3),
            type=int,
            help="Primitive n-th root of unity used to generate GF(65537)",
        )

        self.arg_circ()

        subparsers = self.main_parser.add_subparsers(dest="command")

        backup = subparsers.add_parser(
            "backup",
            parents=[self.common_parser],
            help="Full backup + optional extra parity",
            formatter_class=self.HelpFormatter,
        )
        backup.set_defaults(ecc_ratio=CLI.DEFAULT("1/1", parse_ratio_list))
        backup.add_argument("input_file", metavar="file", nargs="+")

        existing = backup.add_mutually_exclusive_group()
        existing.add_argument("-f", "--force", action="store_true", help="Overwrite existing parity files")
        existing.add_argument("-i", "--ignore-existing", action="store_true", help="Ignore existing parity files")
        existing.add_argument(
            "-u",
            "--update",
            action="store_true",
            help="Update parity files if changes are detected (see --update-check)",
        )
        existing.add_argument(
            "--update-check",
            choices=["timestamp", "hash"],
            default=CLI.DEFAULT("timestamp"),
            help="Update parity files based on timestamp or hash",
        )

    def arg_simd(self):
        simd_choices = ["x4", "x8", "x16", "sse", "avx2", "avx512"]
        self.check(self.check_mutually_exclusive, "--simd", "--no-simd")
        self.common_parser.add_argument(
            "--simd",
            metavar="simd",
            nargs="?",
            const="auto",
            default=CLI.DEFAULT("auto"),
            type=parse_choice_list(["auto"] + simd_choices),
            help="Enable SIMD optimizations",
        )
        self.common_parser.add_argument(
            "--no-simd",
            metavar="simd",
            nargs="?",
            const="*",
            type=parse_choice_list(["*"] + simd_choices),
            help="Disable SIMD optimizations",
        )

    def arg_circ(self):
        circ = self.common_parser.add_argument_group("circ16")

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
            default=CLI.DEFAULT("508,512,1020,1024,2040,2048,4088,4096", parse_size_list),
            type=parse_size_list,
            help="Payload size used by the inner codec",
        )
        add_circ_arg(
            "--inner-ecc",
            metavar="size",
            action="store",
            default=CLI.DEFAULT("4,8", parse_size_list),
            type=parse_size_list,
            help="ECC size of the inner codec",
        )
        add_circ_arg(
            "--inner-interleaved-ecc",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Total ECC size of the inner codec (inner_ecc * outer_block * outer_interleave)",
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
            default=CLI.DEFAULT("16,32,64,128,256", parse_size_list),
            type=parse_size_list,
            help="ECC size of the outer codec",
        )
        add_circ_arg(
            "--outer-interleaved-ecc",
            metavar="size",
            action="store",
            type=parse_size_list,
            help="Total ECC size of the outer codec (outer_ecc * inner_message * outer_interleave)",
        )

        add_circ_arg(
            "--outer-interleave",
            metavar="count",
            action="store",
            # default=CLI.DEFAULT("1", parse_int_list),
            type=parse_int_list,
            help="Interleave multiple blocks",
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
                    setattr(namespace, arg, CLI.Value(None))

    def run_checks(self, namespace):
        for func, args in self.checks:
            func(namespace, *args)

    def parse_args(self):
        args = self.main_parser.parse_args()

        for key, value in args.__dict__.items():
            if not isinstance(value, self.Value):
                setattr(args, key, self.Value(value))

        try:
            self.run_checks(args)
        except ValueError as e:
            # import traceback; traceback.print_exc()
            self.main_parser.error(str(e))

        return args
