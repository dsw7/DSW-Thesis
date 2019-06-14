"""
dsw7@sfu.ca
The "main" script that accepts PDB code input from user.
"""

import os
from sys import path; path.append(r"utils")
from ma import MetAromatic
from pprint import pprint
from argparse import ArgumentParser
from pymongo import MongoClient

COLUMNS = ["ARO", "ARO POS", "MET", "MET POS", "NORM", "MET-THETA", "MET-PHI"]
DEFAULT_PORT = 27017
DEFAULT_HOST = "localhost"
DB = "met_aromatic"
COL = "results"

parser = ArgumentParser(description="Input a PDB code:")
parser.add_argument('--code', help='Usage: $ python runner.py --code <1abc>', default='0', type=str)
parser.add_argument('--path', help='Usage: $ python runner.py --path /path/to/txt', default='0', type=str)
parser.add_argument('--cutoff', help='Usage: $ python runner.py --cutoff <float>', default=6.0, type=float)
parser.add_argument('--angle', help='Usage: $ python runner.py --angle <float>', default=109.5, type=float)
parser.add_argument('--model', help='Usage: $ python runner.py --model <cp|rm>', default='cp')
parser.add_argument('--verbose', help='Usage: $ python runner.py --verbose <cp|rm>', action='store_true')
parser.add_argument('--export', help='Usage: $ python runner.py --export', action='store_true')

code = parser.parse_args().code
path = parser.parse_args().path
cutoff = parser.parse_args().cutoff
angle = parser.parse_args().angle
model = parser.parse_args().model
verbose = parser.parse_args().verbose
export = parser.parse_args().export  # export to MongoDB

client = MongoClient(DEFAULT_HOST, DEFAULT_PORT)
db = client[DB]
col = db[COL]


def verify_user_input():
    if (code == '0') and (path == '0'):
        exit('No PDB code or path to .txt given.')
    elif (len(code) != 4) and (code != '0'):
        exit('Invalid pdb code: {}'.format(code))
    elif (code != '0') and (path != '0'):
        exit('Cannot choose between .txt file and pdb code.')
    else:
        pass


def print_args():
    print("Analyzing: {}".format(code))
    print("Cutoff: {}".format(cutoff))
    print("Angle: {}".format(angle))
    print("Model: {}".format(model))
    print("Mongo Port: {}".format(DEFAULT_PORT))
    print("Mongo Host: {}".format(DEFAULT_HOST))
    print("Database Name: {}".format(DB))
    print("Collection Name: {}\n".format(COL))


def read_pdb_code_txt_file(filepath):
    with open(filepath, 'r') as f:
        codes = f.read().splitlines()

    if not codes:
        exit('Empty .txt file.')
    else:
        return codes


def run_met_aromatic(pdbcode):
    try:
        print('-' * 10)
        print('Code: {}'.format(pdbcode))
        return MetAromatic(pdbcode, cutoff=cutoff, angle=angle, model=model).met_aromatic()
    except Exception as exception:
        print('An exception has occurred:')
        print(exception)


def mapper(result):
    outgoing = []
    for item in result:
        outgoing.append({
            "aro": item[0],
            "arores": item[1],
            "met": item[3],
            "norm": item[4],
            "met-theta": item[5],
            "met-phi": item[6]
        })
    return outgoing


#TODO: continue here with export to mongodb
if __name__ == '__main__':
    verify_user_input()
    print_args()

    if (code != '0') and (path == '0'):  # i.e. user inputs a valid pdb code but no path
        results = run_met_aromatic(code)
        if not results:
            print('No interactions.')
        else:
            results = mapper(results)
            if verbose:
                pprint(results)

    elif (code == '0') and (path != '0'):
        if not os.path.exists(path):
            exit('Path to file does not exist.')
        else:
            pdb_codes = read_pdb_code_txt_file(path)

        for code in pdb_codes:
            results = run_met_aromatic(code)
            if not results:
                print('No interactions.')
            else:
                results = mapper(results)
                if verbose:
                    pprint(results)

client.close()
