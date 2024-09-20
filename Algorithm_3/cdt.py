# Author: Ben Lobo (lobo@virginia.edu)

import os
import shutil

base_path = ''

        
def create_clean_directory(output_path):
    """
    Attempts to create the output_path directory, raises an error
    if the error message is anything other than the directory already
    existing.  Cleans the directory if it already exists.
    @type output_path: str
    """
    create_dir(output_path)  # Create the directory if it does not already exist

    # Clean the directory if it is not already empty
    for output in os.listdir(output_path):
        the_path = f'{output_path}/{output}'
        if os.path.isdir(the_path):
            shutil.rmtree(the_path)
        elif os.path.isfile(the_path):
            os.remove(the_path)
        else:
            assert False, 'Output in folder is neither a file nor a folder'


def create_dir(dir_path):
    """Create the directory at dir_path if it does not exist.
    @type dir_path: str
    """
    try:
        os.makedirs(dir_path)
    except OSError:
        if not os.path.isdir(dir_path):
            raise     