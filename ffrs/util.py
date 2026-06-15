#  util.py
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

import ffrs

_B64HEX_ALPHABET = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-_"


class LazyString:
    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        return str(self.func(*self.args, **self.kwargs))


def b64hex(val: int) -> str:
    result = ""
    while val:
        result = _B64HEX_ALPHABET[val & 0x3F] + result
        val >>= 6
    return result


def b64hex_dec(val: str) -> int:
    result = 0
    for c in val:
        result <<= 6
        result |= _B64HEX_ALPHABET.index(c)
    return result


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
    return s + " " * max(0, w - len(re.sub(r"\033\[[0-9;]*m", "", s)))


def format_config(config) -> str:
    parts = ["config: {"]
    for k in sorted(config):
        val = config[k]
        if val is None:
            sval = "\033[35m[[[...]]]\033[0m"
        elif len(val) > 32:
            sval = f"\033[35m[[[{len(val)}]]]\033[0m"
        else:
            sval = "[" + ", ".join((map(format_size, sorted(val)))) + "]"

        parts.append(f"  {repr(k):23}: {sval},")
    parts.append("}\n")
    return "\n".join(parts)


def format_config_args(config) -> str:
    assert (
        len(config["inner_block"])
        == len(config["inner_ecc"])
        == len(config["outer_block"])
        == len(config["outer_ecc"])
        == len(config["interleave"])
        == 1
    )
    return (
        f"CIRC16({next(iter(config['inner_block'])) // 2}, {next(iter(config['inner_ecc'])) // 2},"
        f" {next(iter(config['outer_block']))}, {next(iter(config['outer_ecc']))}, {next(iter(config['interleave']))})"
    )


def instantiate_config(args, config):
    return ffrs.CIRC16(
        next(iter(config["inner_block"])) // 2,
        next(iter(config["inner_ecc"])) // 2,
        next(iter(config["outer_block"])),
        next(iter(config["outer_ecc"])),
        next(iter(config["interleave"])),
        primitive=args.primitive.get(),
        **args.simd.get(),
    )
