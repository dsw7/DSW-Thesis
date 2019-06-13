"""
dsw7@sfu.ca
The "main" script that accepts PDB code input from user.
"""

from sys import path
path.append(r"utils")
from ma import MetAromatic
from pandas import DataFrame
from argparse import ArgumentParser
COLUMNS = ["ARO", "ARO POS", "MET", "MET POS", "NORM", "MET-THETA", "MET-PHI"]

parser = ArgumentParser(description="Input a PDB code:")
parser.add_argument('--code', help='Usage: $ python runner.py --code <1abc>', default='1rcy')
parser.add_argument('--cutoff', help='Usage: $ python runner.py --cutoff <float>', default=6.0, type=float)
parser.add_argument('--angle', help='Usage: $ python runner.py --angle <float>', default=109.5, type=float)
parser.add_argument('--model', help='Usage: $ python runner.py --model <cp|rm>', default='cp')

code = parser.parse_args().code
cutoff = parser.parse_args().cutoff
angle = parser.parse_args().angle
model = parser.parse_args().model

if len(code) != 4:
    exit('Invalid pdb code: {}'.format(code))
else:
    pass

try:
    print("Analyzing: {}".format(code))
    print("Cutoff: {}".format(cutoff))
    print("Angle: {}".format(angle))
    print("Model: {}\n".format(model))
    results = MetAromatic(code, cutoff=cutoff, angle=angle, model=model).met_aromatic()
    if not results:
        print("No interactions.")
    else:
        df = DataFrame(results)
        df.columns = COLUMNS
        print(df)
except Exception as exception:
    print("An exception has occurred:")
    print(exception)
else:
    print('-' * 65)
