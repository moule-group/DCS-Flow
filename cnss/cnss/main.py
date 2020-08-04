from ase.cli.main import main as ase_main

version = '0.0.1'

commands = [
    ('relax', 'cnss.relax'),
    ('phonons', 'cnss.phonons'),
    ('oclimax', 'cnss.oclimax'),
    ('workflow', 'cnss.workflow')]

def main():
    ase_main('cnss', 'CNSS command-line tool', version=version, commands=commands)
