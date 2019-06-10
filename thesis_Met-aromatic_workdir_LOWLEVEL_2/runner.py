"""
dsw7@sfu.ca
"""

from sys import path
path.append(r"utils")
from ma import MetAromatic
from pandas import read_csv, DataFrame

START = 1500
END = 1600
COLUMNS = ["ARO", "ARO POS", "MET", "MET POS", "NORM", "MET-THETA", "MET-PHI"]

df = read_csv("randomized_pdb_codes.csv")
df = df.iloc[START:END]
df.columns = ["CODE"]

for code in df.CODE:
    try:
        print("Analyzing: {}".format(code))
        results = MetAromatic(code).met_aromatic()
        if not results:
            print("No interactions.")
        else:
            df = DataFrame(results)
            df.columns = COLUMNS
            print(df)
    except Exception as exception:
        print("An exception has occurred:")
        print(exception)
