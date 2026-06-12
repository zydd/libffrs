#  par/fs.py
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

import fnmatch
import os
import pathlib
import re

import ffrs
import ffrs.util

from . import log as parent_log

log = parent_log.getChild("fs")


def translate_rule(rule):
    if rule == r"*":
        return r"[^/]+"
    rule_tr = fnmatch.translate(rule)[:-2]
    rule_tr = re.sub(r"(?<!\\)((?:\\\\)*)\.", r"\1[^/]", rule_tr)
    return rule_tr


def parse_exclude_rules(root_dir, rules, regex=None):
    root_dir = root_dir.replace(".", r"\.").rstrip("/")
    combined_rule = regex.pattern if regex else ""
    for rule in rules:
        rule = rule.rstrip("\n")

        if negated := rule.startswith("!"):
            rule = rule[1:]

        if root := rule.startswith("/"):
            rule = rule[1:]

        if not rule:
            continue

        re_rule = root_dir if root else ".*"
        for sub in rule.split("/"):
            if sub:
                re_rule += "/" + translate_rule(sub)
        re_rule += r"/"

        if not rule.endswith("/"):
            re_rule += r"?"
        re_rule += r"$"

        if negated:
            combined_rule = rf"(?!{re_rule}){combined_rule}"
        else:
            combined_rule = rf"(?:{re_rule})|{combined_rule}" if combined_rule else re_rule
        # print(f"{"!" if negated else " "}{rule:30}", re_rule)

    return re.compile(combined_rule)


def sorted_walk(dir, exclude_rules_file, exclude_rules=re.compile("")):
    root, dirs, files = next(os.walk(dir))
    dirs.sort()
    files.sort()

    if exclude_rules_file in files:
        with open(exclude_rules_file, "r") as f:
            exclude_rules = parse_exclude_rules(root, f.readlines(), exclude_rules)

    for file in files:
        path = os.path.join(root, file)
        if not exclude_rules.fullmatch(path):
            yield path

    for dir in dirs:
        path = os.path.join(root, dir, "")
        if not exclude_rules.fullmatch(path):
            yield from sorted_walk(path, exclude_rules_file, exclude_rules)


def input_files_iter(input_path, exclude_rules_file, output_file):
    output_path = pathlib.Path(output_file).resolve()

    for pathname in input_path:
        path = pathlib.Path(pathname)
        if not path.exists():
            log.critical("input path '%s' does not exist", pathname)
            return 1

        if path.is_dir():
            path = path.resolve()
            exclude_rule = re.compile("")
            if output_path.is_relative_to(path):
                relative_output_path = output_path.relative_to(path)
                exclude_rule = re.compile(rf"^./{str(relative_output_path)}$")

            path_trim = len(pathname) + 1
            for file in sorted_walk(pathname, exclude_rules_file, exclude_rules=exclude_rule):
                yield file[path_trim:]

        elif path.is_file():
            yield pathname

        else:
            raise NotImplementedError("unsupported file type: %s", pathname)


def chunk_filelist(input_files, chunk_size):
    files = []
    total_size = 0
    for filename in input_files:
        stat = os.stat(filename)

        if stat.st_size == 0:
            continue

        total_size += stat.st_size
        files.append(("f", stat.st_size, stat.st_mtime_ns, None, filename))

        offset = 0
        written_size = stat.st_size - (total_size - chunk_size)
        while total_size > chunk_size:
            files[-1] = ("c", written_size, offset, filename)
            yield chunk_size, files

            total_size -= chunk_size
            offset += written_size
            written_size = min(total_size, chunk_size)
            files = [("c", written_size, offset, filename)]

        if files[-1][0] == "c":
            files[-1] = ("fc", stat.st_size, offset, stat.st_mtime_ns, None, filename)

        if total_size == chunk_size:
            yield total_size, files
            files = []
    yield total_size, files


def read_filelist_into(filelist, buffer):
    buf_offset = 0
    for type_, *file_info in filelist:
        if type_ == "f":
            # if fd:
            #     log.debug("close file '%s'", fd.name)
            #     del fd

            file_size, mtime_ns, _hash, path = file_info
            size = file_size
            offset = 0
            # log.debug("read file size %s '%s'", ffrs.util.format_size(size), path)
            # fd = open(path, "rb")
            # stat = os.fstat(fd.fileno())
            # assert stat.st_mtime_ns == mtime_ns, "timestamp mismatch for file '%s'" % path
            # assert stat.st_size == size, "file size mismatch for file '%s'" % path

        elif type_ in ["c", "fc"]:
            if type_ == "c":
                size, offset, path = file_info
            else:
                file_size, offset, mtime_ns, _hash, path = file_info
                size = file_size - offset
                # stat = os.stat(path)
                # assert stat.st_size == file_size, "file size mismatch for file '%s'" % path
                # assert stat.st_mtime_ns >= mtime_ns, "timestamp mismatch for file '%s'" % path

            # log.debug("read chunk %s %s '%s'", ffrs.util.format_size(size), ffrs.util.format_size(offset), path)
            # if offset == 0:
            #     if fd:
            #         log.debug("close file '%s'", fd.name)
            #         del fd
            #     fd = open(path, "rb")
            # else:
            #     fd = open(path, "rb")
            #     fd.seek(offset)
            #     assert fd
            #     assert fd.name == path, (fd.name, path)
            #     assert fd.tell() == offset, (fd.tell(), offset)

        else:
            raise NotImplementedError

        with open(path, "rb") as fd:
            stat = os.fstat(fd.fileno())
            assert stat.st_size == file_size, "file size mismatch for file '%s'" % path
            assert stat.st_mtime_ns >= mtime_ns, "timestamp mismatch for file '%s'" % path

            if offset:
                fd.seek(offset)

            read_size = size
            assert size <= len(buffer) - buf_offset
            while read_size > 0:
                read = fd.readinto(buffer[buf_offset : buf_offset + read_size])
                assert read_size == read
                # TODO: simplify
                if read == 0:
                    break
                buf_offset += read
                read_size -= read

    assert buf_offset <= len(buffer), (buf_offset, len(buffer))
    # return fd


def write_filelist_into(filelist, buffer):
    buf_offset = 0
    for type_, *file_info in filelist:
        if type_ == "f":
            size, mtime_ns, _hash, path = file_info
            offset = 0

        elif type_ in ["c", "fc"]:
            if type_ == "c":
                size, offset, path = file_info
            else:
                file_size, offset, mtime_ns, _hash, path = file_info
                size = file_size - offset

        else:
            raise NotImplementedError

        with open(path, "r+b") as fd:
            if offset:
                fd.seek(offset)

            assert size <= len(buffer) - buf_offset
            written = fd.write(buffer[buf_offset : buf_offset + size])
            assert size == written
            buf_offset += written

            atime_ns = os.fstat(fd.fileno()).st_atime_ns

        os.utime(path, ns=(atime_ns, mtime_ns))

    assert buf_offset <= len(buffer), (buf_offset, len(buffer))
    # return fd


def fill_buffer_gen(input_files, buffer_size, buffers):
    file = None

    for i, (combined_size, filelist) in enumerate(chunk_filelist(input_files, buffer_size)):
        buf = buffers[i & 1]
        file = read_filelist_into(filelist, buf.buf, file)
        log.debug("fill buffer %s", buf)
        yield (filelist, buf)

    del file
