import os
import errno
import sys
from contextlib import contextmanager
from pathlib import Path
from ase.io import jsonio

def mkdir(folder):
    try:
        os.mkdir(folder)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

@contextmanager
def mktempdir():
    import tempfile
    import shutil
    
    temp_dir = tempfile.mkdtemp()
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir)
    
@contextmanager
def chdir(folder):
    dir = os.getcwd()
    os.chdir(str(folder))
    yield
    os.chdir(dir)

@contextmanager    
def out(task):
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
        
    return arg

def write_json(filename, data):
    Path(filename).write_text(jsonio.MyEncoder(indent=4).encode(data))

def read_json(filename):
    dct = jsonio.decode(Path(filename).read_text())
    return dct

def done(mode):
    with open(mode + '.done', 'w') as f:
        f.close()

def isdone(mode):
    if os.path.exists(mode + '.done'):
        return True
    else:
        return False

