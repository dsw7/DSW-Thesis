from sys import path; path.append("./utils")
from ma import MetAromatic

if __name__ == "__main__":
    ma_object = MetAromatic("1rcy")
    results = ma_object.met_aromatic()
    df = DataFrame(results)
    print(df)
