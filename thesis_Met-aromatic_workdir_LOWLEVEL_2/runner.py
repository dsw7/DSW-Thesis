from sys import path; path.append(r"utils")
from ma import MetAromatic
from time import sleep
from pandas import DataFrame, read_csv

columns = ["ARO", "ARO POS", "MET", "MET POS", "NORM", "MET-THETA", "MET-PHI"]

df = read_csv("randomized_pdb_codes.csv")
df = df.head(100)
df.columns = ["CODE"]

if __name__ == "__main__":
    for test_entry in df.CODE:
        try:
            print("-" * 50)
            print(test_entry + "\n")
            results = MetAromatic(test_entry).met_aromatic()
            if results == []
                print("No results.")
                pass
            else:
                df = DataFrame(results)
                df.columns = columns
                print(df)
                sleep(0.25)
        except Exception as exception:
            print(exception)


# Common exceptions:
# [WinError 32] The process cannot access the file because it is being used by another process:
# 'C:\\Users\\David\\Dropbox\\git\\DSW-Thesis\\thesis_Met-aromatic_workdir_LOWLEVEL_2\\pdb1acw.ent'
# 1acw, 2a0l
# Length mismatch: Expected axis has 0 elements, new values have 7 elements 1jk8, 5hso, 3jq6
# 4qz6 could not convert string to float: '75.766-106.225'