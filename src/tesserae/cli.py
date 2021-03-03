def main(argv):
    import importlib  # pylint: disable=C0415
    import argparse  # pylint: disable=C0415
    from . import __version__  # pylint: disable=C0415

    prog_name = "tesserae"

    subcommands = {
        "align": "tesserae.command.align.main",
        "profalign": "tesserae.command.profalign.main",
    }

    parser = argparse.ArgumentParser(
        prog=prog_name,
        description="Graph-based mosaic read alignment and exploration algorithms.",
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s version {__version__}"
    )

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument(
        "-q", "--quiet", help="silence logging except errors", action="store_true"
    )
    verbosity_group.add_argument(
        "-v", "--verbose", help="increase output verbosity", action="store_true"
    )
    verbosity_group.add_argument(
        "-vv", "--veryverbose", help="maximal output verbosity", action="store_true"
    )

    parser.add_argument(
        "subcommand",
        choices=sorted(subcommands.keys()),
        help=f"{prog_name} sub-command",
    )
    parser.add_argument("args", nargs=argparse.REMAINDER, help="sub-command arguments")
    args = parser.parse_args(argv[1:])

    from . import log  # pylint: disable=C0415

    log.configure_logging(args)

    package_string, method_string = subcommands[args.subcommand].rsplit(".", 1)
    module = importlib.import_module(package_string)
    return getattr(module, method_string)(args.args)
