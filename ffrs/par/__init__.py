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


import re
import sys


def format_size(n):
    units = ["", "k", "m", "g", "t", "p"]
    colors = [
        "\033[36m",  # cyan
        "\033[32m",  # green
        "\033[33m",  # yellow
        "\033[31m",  # red
        "\033[35m",  # magenta
        "\033[94m",  # light blue
    ]
    if n <= 0:
        return str(n)
    parts = []
    for eng in range(len(units) - 1, -1, -1):
        shift = eng * 10
        value = 1 << shift
        if n >= value:
            q, n = divmod(n, value)
            if q:
                parts.append(f"{colors[eng]}{q}{units[eng]}\033[0m")
    return "".join(parts)


def _ansi_ljust(s, w):
    return s + " " * max(0, w - len(re.sub(r"\x1b\[[0-9;]*m", "", s)))


def print_config(config):
    print()

    if (
        len(config["inner_block"])
        == len(config["inner_ecc"])
        == len(config["outer_block"])
        == len(config["outer_ecc"])
        == 1
    ):
        print(
            f"# CIRC16({(*config['inner_block'],)[0]}, {(*config['inner_ecc'],)[0]},"
            f" {(*config['outer_block'],)[0]}, {(*config['outer_ecc'],)[0]}, {(*config['outer_interleave'],)[0]})"
        )

    print("config: {")
    for k in sorted(config):
        val = config[k]
        if val is None:
            val = "\033[35m[[[...]]]\033[0m"
        elif len(val) > 32:
            val = f"\033[35m[[[{len(val)}]]]\033[0m"
        else:
            val = "[" + ", ".join((map(format_size, sorted(val)))) + "]"

        print(f"  {repr(k):22}: {val},")
    print("}\n")


def print_warning(message):
    print(f"\x1b[33mWarning: {message}\x1b[0m", file=sys.stderr)


constraints = [
    "2 <= {block}",
    "2 <= {message}",
    "2 <= {ecc}",
    #
    "2 <= {inner_block} <= 65536",
    "2 <= {inner_message} <= 65536",
    "2 <= {inner_ecc} <= 32768",
    #
    "2 <= {outer_block} <= 65536",
    "2 <= {outer_message} <= 65536",
    "2 <= {outer_ecc} <= 32768",
    "1 <= {outer_interleave}",
    #
    "{block} == {message} + {ecc}",
    "{inner_block} == {inner_message} + {inner_ecc}",
    "{outer_block} == {outer_message} + {outer_ecc}",
    #
    "{block} == {inner_block} * {outer_block} * {outer_interleave}",
    "{message} == {inner_message} * {outer_message} * {outer_interleave}",
    "{interleaved_ecc} == {inner_message} * {outer_ecc} * {outer_interleave}",
    "{ecc} == {interleaved_ecc} + {outer_block} * {inner_ecc} * {outer_interleave}",
    # "{ecc} == ({inner_message} * {outer_ecc} + {outer_block} * {inner_ecc} * {outer_interleave})",
    #
    "{outer_message} * {outer_ecc_ratio_num} == {outer_ecc} * {outer_ecc_ratio_den}",
    "{message} * {outer_ecc_ratio_num} == {interleaved_ecc} * {outer_ecc_ratio_den}",
    "{outer_block} * {outer_ecc_ratio_num} == {outer_ecc} * ({outer_ecc_ratio_den} + 1)",
    #
    "{outer_ecc} & ({outer_ecc} - 1) == 0",
    #
    "{message} % {inner_message} == 0",
    "{message} % {outer_message} == 0",
    "{message} % {outer_interleave} == 0",
    "{message} % ({inner_message} * {outer_message}) == 0",
    #
    # "{message} % {ecc} == 0",  # not applicable for CIRC
    "{message} % {interleaved_ecc} == 0",
    "{inner_message} % {inner_ecc} == 0",
    "{outer_message} % {outer_ecc} == 0",
    #
    "{block} % {inner_block} == 0",
    "{block} % {outer_block} == 0",
    "{block} % {outer_interleave} == 0",
    "{block} % ({inner_block} * {outer_block}) == 0",
    #
    "{outer_interleave} == {message} // ({inner_message} * {outer_message})",
    "{outer_interleave} == {block} // ({inner_block} * {outer_block})",
    #
    "{ecc} <= {message}",
    "{inner_ecc} <= {inner_message}",
    "{outer_ecc} <= {outer_message}",
    # "math.gcd({outer_ecc_ratio_num}, {outer_ecc_ratio_den}) == 1",
    #
    "{interleaved_ecc} < {ecc}",
]
