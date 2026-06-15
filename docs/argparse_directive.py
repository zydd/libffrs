from docutils import nodes
from docutils.parsers.rst import Directive, directives


import importlib

from docutils import nodes
from docutils.parsers.rst import Directive, directives


class ArgparseDirective(Directive):
    """
    .. argparse::
       :module: mod.__main__
       :func: create_parser
       :prog: mod
    """

    has_content = False

    option_spec = {
        "module": directives.unchanged_required,
        "func": directives.unchanged_required,
        "prog": directives.unchanged,
    }

    def run(self):
        module_name = self.options["module"]
        func_name = self.options["func"]

        module = importlib.import_module(module_name)
        parser_factory = getattr(module, func_name)
        parser = parser_factory().parser

        prog = self.options.get("prog", parser.prog)

        self.override_formatter_width(parser)

        result = []

        # Top-level section
        title = nodes.title(text=prog)
        section = nodes.section(ids=[nodes.make_id(prog)])
        section += title

        help_text = parser.format_help()
        section += nodes.literal_block(help_text, help_text, language="console")

        result.append(section)

        # Find argparse subparsers action
        subparsers_action = None
        for action in parser._actions:
            if action.__class__.__name__ == "_SubParsersAction":
                subparsers_action = action
                break

        if subparsers_action is None:
            return result

        # Create one section per subcommand
        for name, subparser in subparsers_action.choices.items():
            self.override_formatter_width(subparser)
            sub_section = nodes.section(ids=[nodes.make_id(name)])
            sub_section += nodes.title(text=f"{prog} {name}")

            help_text = subparser.format_help()
            sub_section += nodes.literal_block(help_text, help_text, language="console")

            result.append(sub_section)

        return result

    def override_formatter_width(self, parser):
        def wrapped_formatter_class(*args, _formatter_class=parser.formatter_class, **kwargs):
            kwargs["width"] = 80
            return _formatter_class(*args, **kwargs)

        parser.formatter_class = wrapped_formatter_class


def setup(app):
    app.add_directive("argparse", ArgparseDirective)

    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
