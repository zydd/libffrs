#  par/repair.py
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

"""Repair files using parity data"""

import base64
import hashlib
import json
import os
import io

import ffrs
import ffrs.util

from . import cli, format, fs
from .cli import log


def parse_hash(buf: bytes):
    algo, value = buf.split(b":", 1)
    assert algo == b"b2"
    return base64.urlsafe_b64decode(value + b"==").hex()


def parse_file_info(buf: bytes):
    # b"f s:bsr0 t:1yQIW4R8BWI h:b2:KfqN0Swd05U2VAEQoUcUbqHE1r8 p:test\n"

    params = {}
    for part in buf.split():
        key, value = part.split(b":", 1)
        params[key] = value

    path = params[b"p"].decode("utf-8")
    if path.startswith('"'):
        assert len(path) >= 2 and path.endswith('"'), "invalid path encoding"
    else:
        path = f'"{path}"'

    return {
        "hash": parse_hash(params[b"h"]),
        "path": json.loads(path),
        "size": ffrs.util.b64hex_dec(params[b"s"].decode("ascii")),
        "timestamp": ffrs.util.b64hex_dec(params[b"t"].decode("ascii")),
    }


def parse_ffrs_file(f: io.BufferedReader):
    assert f.readline(7) == b"\x89ffrs\r\n"

    codec = None
    file_info = None
    codecs = {b"circ16": ffrs.CIRC16, b"rsi16": ffrs.RSi16}

    while True:
        line = f.readline(4096)
        assert line[-1:] == b"\n", line
        if len(line) == 1:
            break

        key, args = line.split(b" ", 1)

        if key in codecs:
            assert codec is None, "codec redefinition"
            codec = codecs[key].deserialize(line)
        elif key == b"f":
            assert file_info is None, "file info redefinition"
            file_info = parse_file_info(args)

    return codec, file_info


def main(args):
    filename = args.input_file.get()
    if not os.path.isfile(filename):
        log.critical("input file '%s' does not exist", filename)
        return 1

    with format.Reader(filename) as input:
        rs = input.read_header()
        ecc = ffrs.create_buffer(rs.ecc_size)
        data = ffrs.create_buffer(rs.message_size)

        for filelist in input.blocks(ecc):
            fs.read_filelist_into(filelist, data)
            rs.repair(data, ecc)

            fs.write_filelist_into(filelist, data)

    log.info("done")
    return 0


def arg_parser(parser):
    cli_parser = cli.CLI(parser, main)
    cli_parser.parser.add_argument("input_file", metavar="file")
    return cli_parser


if __name__ == "__main__":
    cli.main()
