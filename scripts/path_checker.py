import os
import warnings
import sys
from os.path import exists, abspath


def path_checker(name, specific_path, general_path):
    """Class to check for the invalid path"""

    print(f"Establishing the path for {name}")
    if specific_path:
        path = specific_path
    elif general_path:
        path = general_path
    else:
        warnings.warn("Directory not specified for '" + name + "' in config via local:path_* or *:*:local_path!. Setting to current directory")
        path = "."

    if not exists(path):
        raise FileNotFoundError("Directory not found for '" + name + "' at: '" + path + "' specified via '*:*:local_path' or 'local:path_*'!")
        sys.exit(1)

    path.rstrip(os.sep)
    path = abspath(path)

    return(path)