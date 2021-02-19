import os
import errno
import sys
from contextlib import contextmanager
from pathlib import Path
from ase.io import jsonio

def mkdir(folder):
    """Creates folder in current wd unless OSError occurs.

    Args:
        folder (str): Name of folder to create.
    """
    try:
        os.mkdir(folder)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

@contextmanager
def mktempdir():
    """Creates a temporary directory. 

    Yields:
        str: Temporary directory name.
    """
    import tempfile
    import shutil
    
    temp_dir = tempfile.mkdtemp()
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir)

@contextmanager
def chdir(folder):
    """Changes the working directory to folder if not already current wd. 

    Args:
        folder (str): Name of folder to make wd.
    """
    dir = os.getcwd()
    os.chdir(str(folder))
    yield
    os.chdir(dir)

@contextmanager    
def out(task):
    """Creates .out and .err files for specified task. 

    Args:
        task (str): Name of task. 
    """
    sys.stdout = open(task + '.out', 'a')
    sys.stderr = open(task + '.err', 'a')
    yield    

def get_default_parameters(func):
    """Inspects function func for argument names and default arguments.

    Args:
        func (func): Input function to inspect.

    Returns:
        dict: Keys are argumnent names, values are the default arguments.
    """
    import inspect

    arg_names = inspect.getfullargspec(func).args
    arg_default = inspect.getfullargspec(func).defaults

    arg = {}
    for n, d in zip(arg_names, arg_default):
        arg['{}' .format(n)] = d
        
    return arg

def write_json(filename, data):
    """Creates json file with specified name and data. 

    Args:
        filename (string): Name of json file. Must include .json. 
        data (dict): Dictionary with data. 
    """
    Path(filename).write_text(jsonio.MyEncoder(indent=4).encode(data))

def read_json(filename):
    """Reads specified json file. 

    Args:
        filename (str): Full name of json file.

    Returns:
        dict: Dictionary with data specified by json file. 
    """
    dct = jsonio.decode(Path(filename).read_text())
    return dct

def done(mode):
    """Create .done file for task as f and closes. 

    Args:
        mode (str): Name of task finished. 
    """
    with open(mode + '.done', 'w') as f:
        f.close()

def isdone(mode):
    """Checks for .done file. 

    Args:
        mode (str): Calculator used for task/specified.

    Returns:
        bool: True if mode.done files exists, false if mode.done doesn't exist.
    """
    if os.path.exists(mode + '.done'):
        return True
    else:
        return False

