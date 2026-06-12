#  par/format.py
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
import json
import os
import re
import shlex

import ffrs
from ffrs.util import b64hex, b64hex_dec

from . import log as parent_log

log = parent_log.getChild("format")


class Writer:
    signature = b"\x89ffrs\r\n"

    def __init__(self, path):
        self.path = path
        self._current_dir = ""

    def __enter__(self):
        self.file = open(self.path, "wb")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()

    def write_header(self, codec):
        codec_info = codec.serialize()
        log.debug(Writer.signature)
        log.debug(codec_info)
        self.file.write(Writer.signature)
        self.file.write(codec_info)
        self.file.write(b"\n")

    def _serialize_path(self, path):
        dir, file = os.path.split(path)

        dir = json.dumps(dir, ensure_ascii=False)
        if not re.findall(r"\s", dir):
            dir = dir[1:-1]

        file = json.dumps(file, ensure_ascii=False)
        if not re.findall(r"\s", file):
            file = file[1:-1]

        if dir != self._current_dir:
            self._current_dir = dir
            log.debug("p p:%s", dir)
            self.file.write(f"p p:{dir}\n".encode("utf8"))

        return file

    def _add_file(self, size, mtime_ns, hash, path):
        filename = self._serialize_path(path)
        log.debug(f"f s:%d m:%d h:b2:%s p:%s", size, mtime_ns, hash.hex(), path)
        self.file.write(f"f s:{b64hex(size)} m:{b64hex(mtime_ns)} h:b2:".encode("ascii"))
        self.file.write(base64.urlsafe_b64encode(hash).strip(b"="))
        self.file.write(f" p:{filename}\n".encode("utf8"))

    def _add_file_chunk(self, size, offset, path):
        filename = self._serialize_path(path)
        log.debug(f"c s:{size} o:{offset} p:{path}")
        self.file.write(f"c s:{b64hex(size)} o:{b64hex(offset)} p:{filename}\n".encode("utf8"))

    def _add_file_chunk_final(self, size, offset, mtime_ns, hash, path):
        filename = self._serialize_path(path)
        log.debug(f"fc s:%d o:%d m:%d h:b2:%s p:%s", size, offset, mtime_ns, hash.hex(), path)
        self.file.write(f"fc s:{b64hex(size)} o:{b64hex(offset)} m:{b64hex(mtime_ns)} h:b2:".encode("ascii"))
        self.file.write(base64.urlsafe_b64encode(hash).strip(b"="))
        self.file.write(f" p:{filename}\n".encode("utf8"))

    def write_block(self, filelist, data):
        for t, *args in filelist:
            if t == "f":
                self._add_file(*args)
            elif t == "c":
                self._add_file_chunk(*args)
            elif t == "fc":
                self._add_file_chunk_final(*args)
            else:
                raise NotImplementedError(t)
        log.debug("d s:%d", len(data))
        self.file.write(f"d s:{b64hex(len(data))}\n".encode("ascii"))
        self.file.write(data)
        self.file.write(b"\n")


class Reader:
    signature = b"\x89ffrs\r\n"

    def __init__(self, path):
        self.path = path
        self.file_size = None
        self._current_dir = ""

    def __enter__(self):
        self.file = open(self.path, "rb")
        self.file_size = os.fstat(self.file.fileno()).st_size
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()

    def read_header(self):
        magic = self.file.readline(16)
        assert magic == Reader.signature
        codec_info = self.file.readline(4096)
        log.debug(codec_info)
        codec = ffrs.CIRC16.deserialize(codec_info)
        return codec

    @staticmethod
    def _deserialize_path(path):
        if path.startswith('"'):
            assert len(path) >= 2 and path.endswith('"'), "invalid path encoding"
        else:
            path = f'"{path}"'
        return json.loads(path)

    @staticmethod
    def _deserialize_hash(buf: bytes):
        algo, value = buf.split(":", 1)
        assert algo == "b2"
        return base64.urlsafe_b64decode(value + "==").hex()

    def _parse_header_line(self, line):
        assert line.endswith("\n")
        type_, args = line.split(" ", 1)

        params = {}
        for arg in shlex.split(args):
            key, value = arg.split(":", 1)
            if key in ["m", "n", "o", "s"]:
                value = b64hex_dec(value)
            elif key == "p":
                value = self._deserialize_path(value)
            elif key == "h":
                value = self._deserialize_hash(value)
            else:
                raise NotImplementedError(key)
            params[key] = value

        return type_, params

    def read_block(self, ecc_buffer):
        filelist = []
        while True:
            line = self.file.readline(4096)
            assert line[-1:] == b"\n", line
            if len(line) == 1:
                break

            key, params = self._parse_header_line(line.decode("utf8"))

            if key == "f":
                filelist.append(
                    ("f", params["s"], params["m"], params["h"], os.path.join(self._current_dir, params["p"]))
                )
            elif key == "c":
                filelist.append(("c", params["s"], params["o"], os.path.join(self._current_dir, params["p"])))
            elif key == "fc":
                filelist.append(
                    (
                        "fc",
                        params["s"],
                        params["o"],
                        params["m"],
                        params["h"],
                        os.path.join(self._current_dir, params["p"]),
                    )
                )
            elif key == "d":
                data_size = params["s"]
                assert data_size == len(ecc_buffer)
                read_size = self.file.readinto(ecc_buffer)
                assert read_size == len(ecc_buffer)
                assert self.file.read(1) == b"\n"
                return filelist
            elif key == "p":
                self._current_dir = params["p"]
            else:
                raise NotImplementedError(key)

    def blocks(self, ecc_buffer):
        while self.file.tell() < self.file_size:
            yield self.read_block(ecc_buffer)
