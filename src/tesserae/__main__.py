"""
Entrypoint module, in case you use `python -mnameless`.
Why does this file exist, and why __main__? For more info, read:
- https://www.python.org/dev/peps/pep-0338/
- https://docs.python.org/2/using/cmdline.html#cmdoption-m
- https://docs.python.org/3/using/cmdline.html#cmdoption-m
"""

import sys

from .cli import main


def main_entry():
    """Main entry point into tesserae"""

    main(sys.argv)


if __name__ == "__main__":
    main_entry()
