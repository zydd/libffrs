#  par/backup.py
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

"""Full backup + optional extra parity"""

import base64
import hashlib
import json
import os
import re

import ffrs
import ffrs.util

from . import cli, format, opt
from .cli import log


def main(args):
    for filename in args.input_file.get():
        if not os.path.isfile(filename):
            log.critical("input file '%s' does not exist", filename)
            return 1

        basename = os.path.basename(filename)
        basename = json.dumps(basename, ensure_ascii=False)
        if not re.findall(r"\s", basename):
            basename = basename[1:-1]
        output_filename = filename + ".ffrs"
        log.info("output file: %s", output_filename)

        if os.path.isfile(output_filename):
            if args.force.get():
                log.info("overwriting '%s'", output_filename)
            else:
                log.critical("output file '%s' already exists", output_filename)
                return 1

        input_stat = os.stat(filename)

        if input_stat.st_size == 0:
            log.debug("skipping empty file: '%s'", filename)
            continue

        rs = opt.maximize_interleaving(args, input_stat.st_size)
        log.info("codec: %s", rs)
        log.info("buffer size: %s", ffrs.util.format_size(rs.message_size))

        # TODO: chunking if file is too large
        assert rs.message_size >= input_stat.st_size

        hasher = hashlib.blake2b(digest_size=20, usedforsecurity=False)
        buffer = ffrs.create_buffer(rs.message_size)
        assert len(buffer) == rs.message_size

        with open(filename, "rb") as input_file:
            while read := input_file.readinto(buffer):
                if read < rs.message_size:
                    buffer[read:] = bytes(rs.message_size - read)

                hasher.update(buffer[:read])
                encoded = rs.encode(buffer)

            with format.Writer(output_filename) as fd:
                fd.write_header(rs)

                file_info = ("f", input_stat.st_size, input_stat.st_mtime_ns, hasher.digest(), filename)
                fd.write_block([file_info], encoded)

    log.info("done")
    return 0


def arg_parser(parser):
    cli_parser = cli.CLI(parser, main)
    cli_parser.parser.set_defaults(ecc_ratio="1/1")
    backup = parser.add_argument_group("backup")
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
        default=cli.DEFAULT("timestamp"),
        help="Update parity files based on timestamp or hash",
    )
    return cli_parser


if __name__ == "__main__":
    cli.main()
