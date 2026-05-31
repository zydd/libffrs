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


import logging
import os
import sys


class FfrsParException(Exception):
    pass


class ColorFormatter(logging.Formatter):
    bold = "\033[1m"
    bold_gray = "\033[38;1m"
    bold_red = "\033[31;1m"
    bold_white = "\033[37;1m"
    bold_yellow = "\033[33;1m"
    dim = "\033[1;30m"
    gray = "\033[38;20m"
    red = "\033[31;20m"
    white = "\033[37;20m"
    yellow = "\033[33;20m"
    reset = "\033[0m"
    # format_string = "%(asctime)s %(threadName)-8s %(levelname)-8s %(filename)s:%(lineno)-3d: %(funcName)s: %(message)s"
    format_string = "[%(levelname)-8s %(relpath)s:%(lineno)d:%(funcName)s] "
    module_path = os.path.dirname(os.path.dirname(__file__))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.FORMATS = {
            logging.DEBUG: logging.Formatter(f"{self.dim}{self.format_string}%(message)s{self.reset}"),
            logging.INFO: logging.Formatter(f"{self.format_string}%(message)s"),
            logging.WARNING: logging.Formatter(
                f"{self.bold_yellow}{self.format_string}{self.reset}{self.yellow}%(message)s{self.reset}"
            ),
            logging.ERROR: logging.Formatter(
                f"{self.bold_red}{self.format_string}{self.reset}{self.red}%(message)s{self.reset}"
            ),
            logging.CRITICAL: logging.Formatter(f"{self.bold_red}{self.format_string}%(message)s{self.reset}"),
        }

    def format(self, record):
        record.relpath = os.path.relpath(record.pathname, self.module_path)
        return self.FORMATS[record.levelno].format(record)


logger = logging.getLogger("par")
logger.setLevel(logging.ERROR)
ch = logging.StreamHandler(stream=sys.stdout)
ch.setFormatter(ColorFormatter())
logger.addHandler(ch)


constraints = [
    "2 <= {block}",
    "2 <= {message}",
    "2 <= {ecc}",
    #
    "2 <= {inner_block} <= 65536",
    "2 <= {inner_message} <= 4096",
    "2 <= {inner_ecc} <= 32768",
    #
    "2 <= {outer_block} <= 65536",
    "2 <= {outer_message} <= 65536",
    "2 <= {outer_ecc} <= 32768",
    "1 <= {interleave}",
    "1 <= {outer_interleave}",
    #
    "127 <= {inner_message} // {inner_ecc} <= 512",
    #
    "{block} == {message} + {ecc}",
    "{inner_block} == {inner_message} + {inner_ecc}",
    "{outer_block} == {outer_message} + {outer_ecc}",
    #
    "{block} == {inner_block} * {outer_block} * {interleave}",
    "{message} == {outer_message} * {outer_interleave}",
    "{outer_interleaved_ecc} == {outer_ecc} * {outer_interleave}",
    "{inner_interleaved_ecc} == {outer_block} * {inner_ecc} * {interleave}",
    "{ecc} == {outer_interleaved_ecc} + {inner_interleaved_ecc}",
    #
    "{outer_message} * {outer_ecc_ratio_num} == {outer_ecc} * {outer_ecc_ratio_den}",
    "{message} * {outer_ecc_ratio_num} == {outer_interleaved_ecc} * {outer_ecc_ratio_den}",
    "{outer_block} * {outer_ecc_ratio_num} == {outer_ecc} * ({outer_ecc_ratio_den} + 1)",
    #
    "{outer_ecc} & ({outer_ecc} - 1) == 0",
    #
    "{message} % {inner_message} == 0",
    "{message} % {outer_message} == 0",
    "{message} % {interleave} == 0",
    "{message} % {outer_interleave} == 0",
    "{message} % ({inner_message} * {outer_message}) == 0",
    # "{message} % {ecc} == 0",  # not applicable for CIRC
    "{message} % {outer_interleaved_ecc} == 0",
    "{inner_message} % {inner_ecc} == 0",
    "{outer_message} % {outer_ecc} == 0",
    #
    "{block} % {inner_block} == 0",
    "{block} % {outer_block} == 0",
    "{block} % {interleave} == 0",
    "{block} % ({inner_block} * {outer_block}) == 0",
    #
    "{interleave} == {message} // ({inner_message} * {outer_message})",
    "{interleave} == {block} // ({inner_block} * {outer_block})",
    "{outer_interleave} == {interleave} * {inner_message}",
]
