import os
import errno
import sys
from contextlib import contextmanager
from pathlib import Path
from ase.io import jsonio

def mkdir(folder):
# make folder, if can't will raise exception --> error and exit
    try:
        os.mkdir(folder)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

@contextmanager
def chdir(folder):
    dir = os.getcwd()
    os.chdir(str(folder))
    yield
    os.chdir(dir)

@contextmanager    
def out(task):
# creates standard out and standard error files (buffer?) that contain results
    sys.stdout = open(task + '.out', 'a')
    sys.stderr = open(task + '.err', 'a')
    yield    

def get_default_parameters(func):
    import inspect
    
    arg_names = inspect.getfullargspec(func).args
    arg_default = inspect.getfullargspec(func).defaults

    arg = {}
    for n, d in zip(arg_names, arg_default):
        arg['{}' .format(n)] = d
        # writing arguments as a dictionary/json
        
    return arg

# writing and reading the json file 
def write_json(filename, data):
    Path(filename).write_text(jsonio.MyEncoder(indent=4).encode(data))

def read_json(filename):
    dct = jsonio.decode(Path(filename).read_text())
    return dct
