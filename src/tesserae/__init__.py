r"""
 _____
|_   _|__  ___ ___  ___ _ __ __ _  ___
  | |/ _ \/ __/ __|/ _ \ '__/ _` |/ _ \
  | |  __/\__ \__ \  __/ | | (_| |  __/
  |_|\___||___/___/\___|_|  \__,_|\___|

             +___________+
            /:\         ,:\
           / : \       , : \
          /  :  \     ,  :  \
         /   :   +-----------+
        +....:../:...+   :  /|
        |\   +./.:...`...+ / |
        | \ ,`/  :   :` ,`/  |
        |  \ /`. :   : ` /`  |
        | , +-----------+  ` |
        |,  |   `+...:,.|...`+
        +...|...,'...+  |   /
         \  |  ,     `  |  /
          \ | ,       ` | /
           \|,         `|/
            +___________+


A pairwise sequence aligner with the ability to jump between
references.
"""

# Set our version:
__version__ = "3.0alpha"

from tesserae.hmm import Tesserae, TesseraeFactory  # type: ignore # noqa: F401
from tesserae.output import SamWriter, pprint_alignment  # noqa: F401
from tesserae.utils import RefInterval, calc_num_mismatches, open_compressed  # noqa: F401
