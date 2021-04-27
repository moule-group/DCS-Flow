from ase.cli.main import main as ase_main

version = '1.0.1'

commands = [
    ('relax', 'dcs.relax'),
    ('md', 'dcs.md'),
    ('chimes', 'dcs.chimes'),
    ('phonons', 'dcs.phonons'),
    ('oclimax', 'dcs.oclimax'),
    ('workflow', 'dcs.workflow'),
    ('train', 'dcs.train')]

def main():
    """Looks up ase_main function, sets up command line.
    """
    ase_main('dcs', 'DCS-Flow command-line tool', version=version, commands=commands)
