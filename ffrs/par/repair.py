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

from . import cli
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
    for filename in args.input_file.get():
        if not os.path.isfile(filename):
            log.critical("input file '%s' does not exist", filename)
            return 1

        # TODO: autodetect file if .ffrs is provided
        assert not filename.endswith(".ffrs")

        input_stat = os.stat(filename)

        if input_stat.st_size == 0:
            log.debug("skipping empty file: '%s'", filename)
            continue

        parity_filename = filename + ".ffrs"
        parity_stat = os.stat(parity_filename)
        with open(parity_filename, "rb") as parity_file:
            rs, file_info = parse_ffrs_file(parity_file)

            log.info("codec: %s", rs)
            log.info("buffer size: %s", ffrs.util.format_size(rs.message_size))
            log.info("file info: %s", file_info)

            assert input_stat.st_mtime_ns == file_info["timestamp"], "timestamp mismatch"
            assert input_stat.st_size == file_info["size"], "file size mismatch"
            assert input_stat.st_size <= rs.message_size, "chunking not supported yet"

            ecc = ffrs.create_buffer(rs.ecc_size)
            read_ecc = parity_file.readinto(ecc)
            assert read_ecc == rs.ecc_size, "could not read ECC data"
            assert parity_file.tell() == parity_stat.st_size, "leftover data in parity file"

        with open(filename, "rb") as input_file:
            data = ffrs.create_buffer(rs.message_size)
            read_data = input_file.readinto(data)

            rs.repair(data, ecc)

            repaired_hash = hashlib.blake2b(data[:read_data], digest_size=20, usedforsecurity=False).hexdigest()
            if repaired_hash != file_info["hash"]:
                log.error("hash mismatch for '%s'", filename)
                return 1

        with open(filename, "wb") as output_file:
            output_file.write(data[:read_data])
        os.utime(filename, ns=(input_stat.st_atime_ns, file_info["timestamp"]))

    log.info("done")
    return 0


def arg_parser(parser):
    cli_parser = cli.CLI(parser, main)
    cli_parser.parser.add_argument("input_file", metavar="file", nargs="+")
    return cli_parser


if __name__ == "__main__":
    cli.main()
