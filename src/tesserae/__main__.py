"""
Entrypoint module, in case you use `python -mnameless`.
Why does this file exist, and why __main__? For more info, read:
- https://www.python.org/dev/peps/pep-0338/
- https://docs.python.org/2/using/cmdline.html#cmdoption-m
- https://docs.python.org/3/using/cmdline.html#cmdoption-m
"""

from __future__ import annotations

import logging
import argparse

import tesserae
from tesserae.cli.registry import SubcommandRegistry
from tesserae.cli.align import AlignSubcommand

logger = logging.getLogger()


class TesseraeCLI(SubcommandRegistry):
    def __init__(self):
        desc = tesserae.__doc__
        desc += f"\n\nVersion: {tesserae.__version__}"

        super().__init__(description=desc, version=tesserae.__version__,
                         formatter_class=argparse.RawDescriptionHelpFormatter)

        # Add arguments to control log level/verboseness
        self.parser.add_argument(
            '-v', '--verbose', action='count', default=0, required=False,
            help="Increase verbosity level."
        )

    def __call__(self, *args, **kwargs):
        args = self.parser.parse_args()

        # Setup logging
        logger.setLevel(logging.WARNING)

        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s:%(name)s:%(message)s")
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        if args.verbose > 0:
            logger.setLevel(logging.INFO)

        if args.verbose > 1:
            logger.setLevel(logging.DEBUG)

        self.run(args)


tesserae_cli = TesseraeCLI()
tesserae_cli.register_subcommand('align', AlignSubcommand())


if __name__ == "__main__":
    tesserae_cli()
