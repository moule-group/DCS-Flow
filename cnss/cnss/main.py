from ase.cli.main import main as ase_main

# define version of cnss?
version = '0.0.2'

# defines commands
commands = [
    ('relax', 'cnss.relax'),
    ('phonons', 'cnss.phonons'),
    ('oclimax', 'cnss.oclimax'),
    ('workflow', 'cnss.workflow')]

def main():
# look up ase_main function, tool that sets up the command line
    ase_main('cnss', 'CNSS command-line tool', version=version, commands=commands)
