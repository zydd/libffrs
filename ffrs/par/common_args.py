#
#  par/common_args.py
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

from . import cli


def add_common_args(cli_inst: cli.CLI):
    add_generic_args(cli_inst)
    add_codec_args(cli_inst)
    add_circ_args(cli_inst)
    add_simd_args(cli_inst)


def add_generic_args(cli_inst: cli.CLI):
    cli_inst.parser.add_argument(
        "-v",
        "--verbosity",
        metavar="level",
        action="append",
        nargs="?",
        const="",
        default=cli.DEFAULT("info", lambda a: [a]),
        type=cli_inst.parse_verbosity,
        help="Increase/set verbosity level",
    )


def add_codec_args(cli_inst: cli.CLI):
    codec_args = cli_inst.parser.add_argument_group("codec")
    codec_args.add_argument(
        "--codec",
        action="store",
        default="circ16",
        choices=["rsi16", "circ16"],
        help="Method use to compute parity",
    )

    cli_inst.check(lambda args: args.codec != "rsi16" or cli_inst.check_either(args, "--block", "--message"))
    cli_inst.check(cli_inst.check_either, "--block", "--message")
    # cli_inst.check(cli_inst.check_mutually_exclusive, "--block", "--message")
    codec_args.add_argument(
        "-b",
        "--block",
        metavar="size",
        action="store",
        type=cli.parse_size_list,
        help="Block size to use when encoding. block = message + ecc",
    )
    codec_args.add_argument(
        "-m",
        "--message",
        metavar="size",
        action="store",
        type=cli.parse_size_list,
        help="Payload size in a parity block",
    )

    cli_inst.check(cli_inst.check_mutually_exclusive, "--ecc", "--ecc-ratio")
    codec_args.add_argument(
        "-e",
        "--ecc",
        metavar="size",
        action="store",
        type=cli.parse_size_list,
        help="ECC size in a parity block",
    )
    codec_args.add_argument(
        "-r",
        "--ecc-ratio",
        metavar="ratio",
        action="store",
        default=cli.DEFAULT("1/15,1/16", cli.parse_ratio_list),
        type=cli.parse_ratio_list,
        help="ECC size as a fraction of the message size",
    )

    codec_args.add_argument(
        "--interleave",
        metavar="size",
        action="store",
        # default=cli.DEFAULT("1", cli.parse_size_list),
        type=cli.parse_size_list,
        help="Interleave multiple blocks",
    )

    codec_args.add_argument(
        "--primitive",
        metavar="int",
        action="store",
        default=cli.DEFAULT(3),
        type=int,
        help="Primitive n-th root of unity used to generate GF(65537)",
    )


SIMD_MAP = {
    "x4": "simd_x4",
    "x8": "simd_x8",
    "x16": "simd_x16",
    "sse": "simd_x4",
    "avx2": "simd_x8",
    "avx512": "simd_x16",
}


def _parse_simd(enable, choices):
    def parser(value):
        if value == "*":
            return {
                "simd_x4": enable,
                "simd_x8": enable,
                "simd_x16": enable,
            }
        elif value == "auto":
            return {
                "simd_x4": None,
                "simd_x8": None,
                "simd_x16": None,
            }
        else:
            res = dict()
            for val in value.split(","):
                res[SIMD_MAP[val]] = enable
            return res

    parser.choices = choices
    return parser


def add_simd_args(cli_inst: cli.CLI):
    simd = cli_inst.parser.add_argument_group("simd")
    # cli_inst.check(cli_inst.check_mutually_exclusive, "--simd", "--no-simd")
    enable_choices = ["auto", "*"] + list(SIMD_MAP.keys())
    disable_choices = ["*"] + list(SIMD_MAP.keys())
    simd.add_argument(
        "--simd",
        metavar="simd",
        action=cli.UpdateAction,
        nargs="?",
        const="auto",
        default=cli.DEFAULT("auto", _parse_simd(True, enable_choices)),
        type=_parse_simd(True, enable_choices),
        help="Enable SIMD optimizations",
    )
    simd.add_argument(
        "--no-simd",
        metavar="simd",
        dest="simd",
        action=cli.UpdateAction,
        nargs="?",
        const="*",
        type=_parse_simd(False, disable_choices),
        help="Disable SIMD optimizations",
    )


def add_circ_args(cli_inst: cli.CLI):
    circ = cli_inst.parser.add_argument_group("circ16")

    def add_circ_arg(*args, **kwargs):
        if "dest" in kwargs:
            dest = kwargs["dest"]
        else:
            assert args[-1].startswith("--")
            dest = args[-1][2:].replace("-", "_")

        cli_inst.check(
            cli_inst.warn,
            lambda args: args.codec.get() == "circ16" or getattr(args, dest).is_default(),
            "circ16 options have no effect when codec != circ16",
        )
        circ.add_argument(*args, **kwargs)

    add_circ_arg(
        "--inner-block",
        metavar="size",
        action="store",
        type=cli.parse_size_list,
        help="Block size used by the inner codec",
    )
    add_circ_arg(
        "--inner-message",
        metavar="size",
        action="store",
        default=cli.DEFAULT("508,512,1020,1024,2040,2048,4088,4096", cli.parse_size_list),
        type=cli.parse_size_list,
        help="Payload size used by the inner codec",
    )
    add_circ_arg(
        "--inner-ecc",
        metavar="size",
        action="store",
        default=cli.DEFAULT("4,8", cli.parse_size_list),
        type=cli.parse_size_list,
        help="ECC size of the inner codec",
    )
    add_circ_arg(
        "--inner-interleaved-ecc",
        metavar="size",
        action="store",
        type=cli.parse_size_list,
        help="Total ECC size of the inner codec (inner_ecc * outer_block * interleave)",
    )

    add_circ_arg(
        "--outer-block",
        metavar="size",
        action="store",
        type=cli.parse_size_list,
        help="Block size used by the outer codec",
    )
    add_circ_arg(
        "--outer-message",
        metavar="size",
        action="store",
        type=cli.parse_size_list,
        help="Payload size used by the outer codec",
    )
    add_circ_arg(
        "--outer-ecc",
        metavar="size",
        action="store",
        default=cli.DEFAULT("16,32,64,128,256", cli.parse_size_list),
        type=cli.parse_size_list,
        help="ECC size of the outer codec",
    )
    add_circ_arg(
        "--outer-interleaved-ecc",
        metavar="size",
        action="store",
        type=cli.parse_size_list,
        help="Total ECC size of the outer codec (outer_ecc * inner_message * interleave)",
    )
    add_circ_arg(
        "--outer-interleave",
        metavar="count",
        action="store",
        # default=cli.DEFAULT("1", cli.parse_int_list),
        type=cli.parse_int_list,
        help="Number of interleaved blocks of the outer codec",
    )
