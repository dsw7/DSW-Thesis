from sys import path; path.append("./utils")
from ma import MetAromatic
from pandas import DataFrame

columns = ["ARO", "ARO POS", "MET", "MET POS", "NORM", "MET-THETA", "MET-PHI"]

if __name__ == "__main__":
    ma_object = MetAromatic("1rcy")
    results = ma_object.met_aromatic()
    df = DataFrame(results)
    df.columns = columns
    print(df)

# TODO: run unit test against grad school version