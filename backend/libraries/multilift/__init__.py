__author__ = "George Young"
__email__ = ""

__prog__ = "multilift"
__version__ = "0.3"
__status__ = "Development"
__website__ = "https://github.com/dbauerlab/multilift"
__prog_string__ = f"{__prog__} v{__version__} ({__status__})"

__license__ = "MIT license"

from .main import get_multilift_sequences, generate_multilift_sequences, multilift

__all__ = ["get_multilift_sequences", "generate_multilift_sequences", "multilift"]
