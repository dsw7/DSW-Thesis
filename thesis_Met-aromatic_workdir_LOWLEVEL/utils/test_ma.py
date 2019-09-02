from ma import MetAromatic
from pandas import DataFrame, read_csv, testing


CUTOFF = 4.9
ANGLE = 109.5
CHAIN = 'A'
COLUMNS = ['ARO', 'ARO RES', 'MET', 'MET RES', 'NORM', 'MET-THETA', 'MET-PHI']
NAMES_CSV = ['ARO', 'ARO RES', 'MET', 'MET RES', 'MET-PHI', 'MET-THETA', 'NORM', 'PDBCODE']
DF_FULL = read_csv('./483OutputA3-3-M-Benchmark.csv', names=NAMES_CSV)


def get_control_dataset(code):
    df_control = DF_FULL[DF_FULL.PDBCODE == code]
    df_control = df_control.drop(['PDBCODE'], axis=1)
    df_control = df_control.reindex(sorted(df_control.columns), axis=1)
    df_control = df_control.sort_values(by='NORM')
    df_control = df_control.reset_index(drop=True)
    return df_control


def get_test_dataset(code):
    ma = MetAromatic(code, chain=CHAIN, cutoff=CUTOFF, angle=ANGLE, model='cp').met_aromatic()
    df_test = DataFrame(ma, columns=COLUMNS)
    df_test = df_test.reindex(sorted(df_test.columns), axis=1)
    df_test = df_test.sort_values(by='NORM')
    df_test = df_test.reset_index(drop=True)
    df_test = df_test.astype({'MET RES': 'int64', 'ARO RES': 'int64'})
    return df_test
    

def test_equality():
    for code in set(DF_FULL.PDBCODE):
        print(code)
        # for code in df_full.PDBCODE.sample(SAMPLE_SIZE):
        df_control = get_control_dataset(code=code)
        df_test = get_test_dataset(code=code)
        print(df_control)
        print(df_test)
        testing.assert_frame_equal(df_test, df_control)
   

if __name__ == '__main__':
    test_equality()

