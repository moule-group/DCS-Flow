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

# def gen2xyz(genfile):
#     import numpy as np
    
#     with open(genfile, 'r') as fp:
#         alllines = fp.readlines()

#     # strip out comments starting with hashmarks
#     lines = []
#     for line in alllines:
#         li = line.partition('#')[0]
#         li = li.rstrip()
#         if li:
#             lines.append(li)

#     words = lines[0].split()
#     natom = int(words[0])
#     flag = words[1].lower()
#     if flag == "s":
#         periodic = True
#     elif flag == "f":
#         periodic = True
#     else:
#         periodic = False
#     specienames = lines[1].split()
#     indexes = np.empty((natom, ), dtype=int)
#     coords = np.empty((natom, 3), dtype=float)
#     for ii, line in enumerate(lines[2:2+natom]):
#         words = line.split()
#         indexes[ii] = int(words[1]) - 1
#         coords[ii] = np.array(words[2:5], dtype=float)
#     if periodic:
#         origin = np.array(lines[natom+2].split(), dtype=float)
#         latvecs = np.empty((3, 3), dtype=float)
#         for jj in range(3):
#             latvecs[jj] = np.array(lines[natom+3+jj].split(), dtype=float)
#     else:
#         origin = None
#         latvecs = None

#     xyzfile = genfile[:-4] + '.xyz'
#     with open(xyzfile, 'w') as f:
#         f.write("{0:d}\n".format(natom))
#         for ii in range(natom):
#             f.write("{0:3s} {1:18.10E} {2:18.10E} {3:18.10E}\n".format(
#             specienames[indexes[ii]], *coords[ii]))
