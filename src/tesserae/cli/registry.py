"""
CLI subcommand registry
"""

from __future__ import annotations

import sys
import textwrap
import argparse
from abc import ABCMeta, abstractmethod
from typing import Optional


class Subcommand(metaclass=ABCMeta):
    """Represents a subcommand with its own argument parser arguments."""

    def register_arguments(self, subparser: argparse.ArgumentParser):
        """This function should register all required arguments for this
        subcommand."""

    @abstractmethod
    def __call__(self, *args, **kwargs):  # type: ignore
        """When the subcommand is used on the command line, this function
        will be called."""


class SubcommandRegistry:
    def __init__(self, version=None, subcommands_title="", parent_parser=None, *args, **kwargs):
        if parent_parser is None:
            self.parser = argparse.ArgumentParser(*args, **kwargs)
            self.parser.set_defaults(subcommand_func=None)
        else:
            self.parser = parent_parser

        self.subparsers = self.parser.add_subparsers(
            title=subcommands_title if subcommands_title else "Subcommands")

        if version:
            self.parser.add_argument(
                '--version', action='version', version=version)

    def register_subcommand(self, name: str, subcommand: Subcommand, **kwargs):
        # Use subcommand class level documentation also for documentation on
        # command line -h/--help
        if hasattr(subcommand.__class__, '__doc__'):
            subcommand_doc = subcommand.__class__.__doc__
            if subcommand_doc:
                first_help_line = subcommand_doc.strip().split('\n\n')[0].strip()
                kwargs['help'] = first_help_line
                kwargs['description'] = textwrap.dedent(subcommand_doc)
                kwargs['formatter_class'] = argparse.RawDescriptionHelpFormatter

        subparser = self.subparsers.add_parser(name, **kwargs)

        # Initialize subcommand arguments
        subcommand.register_arguments(subparser)
        subparser.set_defaults(subcommand_func=subcommand)

    def run(self, parser_args: argparse.Namespace):
        args_dict = vars(parser_args)
        subcommand_func = args_dict.pop('subcommand_func')

        if subcommand_func:
            rc = subcommand_func(**args_dict)
        else:
            self.parser.print_help()
            rc = 1

        if rc is None:
            rc = 0

        sys.exit(rc)


class NestedSubcommandRegistry(Subcommand):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.subparser = None
        self.registry: Optional[SubcommandRegistry] = None

    def register_arguments(self, subparser: argparse.ArgumentParser):
        self.subparser = subparser
        self.registry = SubcommandRegistry(parent_parser=self.subparser)

    def register_subcommand(self, name: str, subcommand: Subcommand, **kwargs):
        if self.registry is None:
            raise ValueError("No registry built yet, use register_arguments first.")

        return self.registry.register_subcommand(name, subcommand, **kwargs)

    def __call__(self, *args, **kwargs):  # type: ignore
        pass
