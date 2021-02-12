from ase.cli.main import main as ase_main

version = '0.0.4'

commands = [
    ('relax', 'cnss.relax'),
    ('md', 'cnss.md'),
    ('chimes', 'cnss.chimes'),
    ('phonons', 'cnss.phonons'),
    ('oclimax', 'cnss.oclimax'),
    ('workflow', 'cnss.workflow'),
    ('train', 'cnss.train')]

def main():
    ase_main('cnss', 'CNSS command-line tool', version=version, commands=commands)
