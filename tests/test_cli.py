#  test_cli.py
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

import hashlib
import logging
import os
import pathlib

import pytest

import ffrs.par.__main__
import ffrs.par.cli
from ffrs.reference.util import randbytes

ffrs.par.cli.log.setLevel(logging.DEBUG)


class TestCli:
    def test_main_help(self, capsys, caplog):
        with pytest.raises(SystemExit) as exc:
            ffrs.par.__main__.main("-h")
        assert exc.value.code == 0
        captured = capsys.readouterr()
        assert captured.out.startswith("usage: ffrs.par")
        assert captured.err == ""
        assert all(record.levelno < logging.WARNING for record in caplog.records)

    @pytest.mark.parametrize("module", ffrs.par.__main__.CLI_MODULES.keys())
    def test_module_help(self, module, capsys, caplog):
        with pytest.raises(SystemExit) as exc:
            ffrs.par.__main__.main(module, "-h")
        assert exc.value.code == 0
        captured = capsys.readouterr()
        assert captured.out.startswith(f"usage: ffrs.par {module}")
        assert captured.err == ""
        assert list(filter(lambda record: record.levelno >= logging.WARNING, caplog.records)) == []

    @pytest.mark.parametrize("input_size", ["128k", "1000000", "1m", "1m256k", "3000k", "3m"])
    def test_benchmark(self, input_size, caplog):
        result = ffrs.par.__main__.main(
            "benchmark", "encode", "--input-size", input_size, "--cooldown=0", "--duration=0"
        )
        assert result == 0
        assert list(filter(lambda record: record.levelno >= logging.WARNING, caplog.records)) == []


@pytest.mark.parametrize("input_size", ["8000", "8k", "64k", "1000000", "1m", "1m256k", "3000k", "3m"])
class TestCliBackupRepair:
    def _backup(self, input_size, tmp_path: pathlib.Path):
        input_size_value = ffrs.par.cli.parse_size(input_size)
        input_file = tmp_path / "input.bin"
        input_file.write_bytes(randbytes(input_size_value))
        assert ffrs.par.__main__.main("backup", str(input_file)) == 0

        output_file = input_file.with_suffix(input_file.suffix + ".ffrs")
        assert output_file.exists()
        max_ratio_diff = 0.15 if input_size_value <= 8192 else 0.05 if input_size_value <= 2**16 else 0.02
        assert abs(output_file.stat().st_size / input_size_value - 1) < max_ratio_diff

        return input_file, output_file

    def test_backup(self, input_size, capsys, caplog, tmp_path: pathlib.Path):
        self._backup(input_size, tmp_path)
        assert capsys.readouterr().err == ""
        assert list(filter(lambda record: record.levelno >= logging.WARNING, caplog.records)) == []

    def test_backup_repair_no_errors(self, input_size, capsys, caplog, tmp_path: pathlib.Path):
        input_file, _output_file = self._backup(input_size, tmp_path)
        input_file_mtime = input_file.stat().st_mtime_ns
        input_file_hash = hashlib.blake2b(input_file.read_bytes()).digest()

        assert ffrs.par.__main__.main("repair", str(input_file)) == 0

        assert input_file.stat().st_mtime_ns == input_file_mtime
        assert hashlib.blake2b(input_file.read_bytes()).digest() == input_file_hash

        assert capsys.readouterr().err == ""
        assert list(filter(lambda record: record.levelno >= logging.WARNING, caplog.records)) == []

    def test_backup_repair_corrupted(self, input_size, capsys, caplog, tmp_path: pathlib.Path):
        input_file, _output_file = self._backup(input_size, tmp_path)
        input_file_mtime = input_file.stat().st_mtime_ns
        input_file_hash = hashlib.blake2b(input_file.read_bytes()).digest()

        # corrupt input file
        input_file.write_bytes(randbytes(input_file.stat().st_size))
        os.utime(input_file, ns=(input_file.stat().st_atime_ns, input_file_mtime))
        assert hashlib.blake2b(input_file.read_bytes()).digest() != input_file_hash

        assert ffrs.par.__main__.main("repair", str(input_file)) == 0

        assert input_file.stat().st_mtime_ns == input_file_mtime
        assert hashlib.blake2b(input_file.read_bytes()).digest() == input_file_hash

        assert capsys.readouterr().err == ""
        assert list(filter(lambda record: record.levelno >= logging.WARNING, caplog.records)) == []
