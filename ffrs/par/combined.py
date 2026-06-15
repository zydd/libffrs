#  par/combined.py
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

"""Combined parity data for multiple files"""

import base64
import concurrent.futures
import hashlib
import json
import multiprocessing.managers
import re

from . import cli, format, fs, opt
from .cli import log


def write_result(args, hash, encoded, buffer):
    log.debug("write results for buffer %s", buffer)
    log.debug("encoded size: %s", len(encoded))
    with open(args.output.get(), "ab") as f:
        for filename, file_hash in hash:
            filename = json.dumps(filename, ensure_ascii=False)
            if not re.findall(r"\s", filename):
                filename = filename[1:-1]
            f.write(b"f h:b2:")
            f.write(base64.urlsafe_b64encode(file_hash).strip(b"="))
            f.write(f" p:{filename}\n".encode("utf-8"))
        f.write(b"\n")


def init_hash_process():
    global hash_class
    global hash_obj
    global hash_file
    hash_class = lambda *a, **kw: hashlib.blake2b(*a, **kw, digest_size=20, usedforsecurity=False)
    hash_obj = None
    hash_file = None


def hash_buffer(buffer, filelist):
    global hash_obj
    global hash_file
    global hash_class

    log.debug("hash buffer %s", buffer)

    buf_offset = 0

    filelist_hash = []

    for type_, *file_info in filelist:
        if type_ == "f":
            assert hash_file is None
            assert hash_obj is None

            size, mtime_ns, _, path = file_info
            hash = hash_class(buffer.buf[buf_offset : buf_offset + size]).digest()
            filelist_hash.append((type_, size, mtime_ns, hash, path))
            buf_offset += size

        elif type_ == "c":
            size, offset, path = file_info

            if offset == 0:
                assert hash_obj is None
                hash_obj = hash_class()
            else:
                assert hash_obj is not None

            hash_obj.update(buffer.buf[buf_offset : buf_offset + size])
            filelist_hash.append((type_, *file_info))
            buf_offset += size

        elif type_ == "fc":
            assert hash_obj is not None

            size, offset, mtime_ns, _, path = file_info

            chunk_size = size - offset
            hash_obj.update(buffer.buf[buf_offset : buf_offset + chunk_size])
            hash = hash_obj.digest()
            hash_obj = None

            filelist_hash.append((type_, size, offset, mtime_ns, hash, path))
            buf_offset += chunk_size

    assert buf_offset <= buffer.size
    return filelist_hash


def init_encode_process(rs):
    global encoder_obj
    encoder_obj = rs
    log.info("codec: %s", encoder_obj)


def encode_buffer(buffer):
    global encoder_obj
    log.debug("encode buffer %s", buffer)
    return encoder_obj.encode(buffer.buf)


def main(args):
    rs = opt.circ(args)

    with (
        concurrent.futures.ProcessPoolExecutor(
            max_workers=1, initializer=init_encode_process, initargs=(rs,)
        ) as encode_process,
        concurrent.futures.ProcessPoolExecutor(max_workers=1, initializer=init_hash_process) as hash_process,
        multiprocessing.managers.SharedMemoryManager() as shm,
        format.Writer(args.output.get()) as output,
    ):
        buf1 = shm.SharedMemory(rs.message_size)
        buf2 = shm.SharedMemory(rs.message_size)

        prev_futures = None

        output.write_header(rs)

        input_files = fs.input_files_iter(args.input_path.get(), args.exclude_rules.get(), args.output.get())
        for filelist, buffer in fs.fill_buffer_gen(input_files, rs.message_size, (buf1, buf2)):
            encode_future = encode_process.submit(encode_buffer, buffer)
            hash_future = hash_process.submit(hash_buffer, buffer, filelist)

            if prev_futures:
                prev_hash_future, prev_encode_future, prev_buf = prev_futures
                output.write_block(prev_hash_future.result(), prev_encode_future.result())

            prev_futures = (hash_future, encode_future, buffer)

        if prev_futures:
            prev_hash_future, prev_encode_future, prev_buf = prev_futures
            output.write_block(prev_hash_future.result(), prev_encode_future.result())

    log.info("done")
    return 0


def create_parser(parent=None):
    parser = cli.create_parser(__loader__.name, __doc__, parent)
    from .common_args import add_common_args

    add_common_args(parser)

    combined = parser.parser.add_argument_group("combined parity")

    combined.add_argument("input_path", metavar="path")
    combined.add_argument("--exclude-rules", metavar="file", default=cli.DEFAULT(".gitignore"))
    combined.add_argument("-o", "--output", metavar="file", help="output file for combined parity data")
    return parser


if __name__ == "__main__":
    cli.main()
