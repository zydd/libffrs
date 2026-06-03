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

import base64
import hashlib
import json
import math
import os
import re
import time

import ffrs
import ffrs.par.exc
import ffrs.util

from . import constraints, logger as parent_logger
from .solver import Solver

logger = parent_logger.getChild("backup")


def maximize_interleaving(args, file_size):
    t0 = time.time()
    file_size = max(file_size, 8192)

    config = dict(
        block=args.block.get(),
        message=None,  # choose based on file_size
        ecc=args.ecc.get(),
        inner_block=args.inner_block.get(range(2, 65536 + 1, 2)),
        inner_message=range(2, 65536 + 1, 2),  # choose based on file_size
        inner_ecc=args.inner_ecc.get(set(2**i for i in range(1, 16))),
        outer_block=args.outer_block.get(range(2, 65536 + 1, 2)),
        outer_message=args.outer_message.get(range(2, 65536 + 1, 2)),
        outer_ecc=args.outer_ecc.get(set(2**i for i in range(1, 16))),
        outer_ecc_ratio_num=(
            set(num for num, den in args.ecc_ratio.get()) if args.ecc_ratio.has_value() else range(1, 1 + 1)
        ),
        outer_ecc_ratio_den=(
            set(den for num, den in args.ecc_ratio.get()) if args.ecc_ratio.has_value() else range(1, 1024 + 1)
        ),
        interleave=args.interleave.get(),
        inner_interleaved_ecc=args.outer_interleaved_ecc.get(),
        outer_interleaved_ecc=args.outer_interleaved_ecc.get(),
        outer_interleave=args.outer_interleave.get(),
    )

    free_variables = {k for k in config if k not in args or not getattr(args, k).has_value()}
    free_variables |= {"message", "inner_message"}
    if args.ecc_ratio.has_value():
        free_variables.remove("outer_ecc_ratio_den")
        free_variables.remove("outer_ecc_ratio_num")

    logger.debug(ffrs.util.format_config(config))

    constraints.append(f"{file_size} <= {{message}} <= {2 * file_size}")

    solver = Solver(config.keys(), constraints, free_variables)

    try:
        solver.solve(config)
    except ValueError:
        logger.exception("could not determine configuration")
        # return
        raise

    for outer_ecc in config["outer_ecc"]:
        try:
            logger.debug("branch outer_ecc")
            config2 = solver.branch_config(config, "outer_ecc", outer_ecc)
            logger.debug(ffrs.util.format_config(config2))

            logger.debug("compute outer_interleave")
            max_padding = 0.05 if file_size <= 16000 else 0.01
            config2["outer_interleave"] = [
                math.ceil(file_size / outer_message) + i
                for outer_message in config2["outer_message"]
                for i in range(min(100, math.ceil(file_size * max_padding / outer_message)))
            ]
            solver.solve(config2, ["outer_interleave"])
            logger.debug(ffrs.util.format_config(config2))

            config2 = solver.maximize(config2, "inner_message")
            logger.debug(ffrs.util.format_config(config2))

            config2 = solver.minimize(config2, "inner_ecc")
            logger.debug(ffrs.util.format_config(config2))

            config2 = solver.minimize(config2, "message")
            logger.debug(ffrs.util.format_config(config2))

            logger.info("solver iterations: %s", solver._constraint_evaluations)
            logger.info("max message ratio: %s", max(config2["message"]) / file_size)
            logger.info("max hash ratio: %s", max(config2["inner_interleaved_ecc"]) / file_size)
            config = config2
            break
        except ValueError:
            logger.debug("could not set outer_ecc to %s", outer_ecc)
    else:
        raise ffrs.par.exc.OptimizationError("could not optimize interleaving")

    logger.debug(ffrs.util.format_config(config))
    t1 = time.time() - t0
    if t1 > 10e-3:
        logger.warning("time to solve: %.2f [ms]", t1 * 1000)
    return ffrs.util.instantiate_config(config)


def main(args):
    for filename in args.input_file.get():
        if not os.path.isfile(filename):
            logger.critical("input file '%s' does not exist", filename)
            return 1

        basename = os.path.basename(filename)
        basename = json.dumps(basename, ensure_ascii=False)
        if not re.findall(r"\s", basename):
            basename = basename[1:-1]
        output_filename = filename + ".ffrs"

        if os.path.isfile(output_filename):
            if args.force.get():
                logger.info("overwriting '%s'", output_filename)
            else:
                logger.critical("output file '%s' already exists", output_filename)
                return 1

        input_stat = os.stat(filename)

        if input_stat.st_size == 0:
            logger.debug("skipping empty file: '%s'", filename)
            continue

        rs = maximize_interleaving(args, input_stat.st_size)
        logger.info("codec: %s", rs)
        logger.info("buffer size: %s", ffrs.util.format_size(rs.message_size))

        # TODO:
        assert rs.message_size >= input_stat.st_size, "Chunking not implemented yet"

        sha1 = hashlib.sha1()
        buffer = ffrs.create_buffer(rs.message_size)
        assert len(buffer) == rs.message_size

        with open(output_filename, "wb") as f:
            f.write(b"\x89ffrs\r\n")
            f.write(rs.serialize())
            f.write(f"\nf s:{ffrs.util.b64hex(input_stat.st_size)}".encode("ascii"))
            f.write(f" t:{ffrs.util.b64hex(input_stat.st_mtime_ns)}".encode("ascii"))
            f.write(b" h:")
            hash_offset = f.tell()
            f.write(base64.urlsafe_b64encode(bytes(sha1.digest_size)).strip(b"="))
            f.write(f" p:{basename}\n\n".encode("utf-8"))

            with open(filename, "rb") as input_file:
                while read := input_file.readinto(buffer):
                    if read < rs.message_size:
                        buffer[read:] = bytes(rs.message_size - read)

                    sha1.update(buffer)
                    encoded = rs.encode(buffer)
                    f.write(encoded)
                f.seek(hash_offset)
                f.write(base64.urlsafe_b64encode(sha1.digest()).strip(b"="))
    logger.info("done")
    return 0
