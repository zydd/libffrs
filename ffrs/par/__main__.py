#  par/__main__.py
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

"""File parity tools"""

import ffrs
import ffrs.par

from . import cli
from .cli import log

import ffrs.par.backup
import ffrs.par.benchmark
import ffrs.par.combined
import ffrs.par.repair

CLI_MODULES = {
    "backup": ffrs.par.backup,
    "benchmark": ffrs.par.benchmark,
    "combined": ffrs.par.combined,
    "repair": ffrs.par.repair,
}


def main(args):
    return args.cli_main.get()(args)


def create_parser(parent=None):
    parser = cli.create_parser(__package__, __doc__, parent)
    subparsers = parser.parser.add_subparsers(dest="command", required=True)
    for name, module in CLI_MODULES.items():
        subparser = subparsers.add_parser(
            name,
            help=module.__doc__.lstrip().split("\n", 1)[0],
            description=module.__doc__,
            formatter_class=parser.parser.formatter_class,
        )
        module.create_parser(subparser)
        subparser.set_defaults(cli_main=module.main)

    return parser


if __name__ == "__main__":
    cli.main()
