"""
Module for decompressing gzip files.
"""

import os
import subprocess

def decompress(files):
    """
    Decompresses a list of gzip files in their current location.

    Args:
        files (list): A list of gzip file paths.

    Returns:
        None
    """

    for f in files:
        print(f)
        try:
            if(os.path.exists(f)):
                directory = os.path.dirname(f)
                subprocess.run(['gzip', '-df', f], check=True, cwd=directory)
            else:
                print(f"File not found for decompression:{f}")
        except subprocess.CalledProcessError as e:
            print(f"Error decompressing {f}: {e}")