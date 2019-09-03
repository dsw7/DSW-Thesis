from ma import MetAromatic
from pandas import DataFrame, read_csv, testing
import pytest


CUTOFF = 4.9
ANGLE = 109.5
CHAIN = 'A'
COLUMNS = ['ARO', 'ARO RES', 'MET', 'MET RES', 'NORM', 'MET-THETA', 'MET-PHI']
NAMES_CSV = ['ARO', 'ARO RES', 'MET', 'MET RES', 'MET-PHI', 'MET-THETA', 'NORM', 'PDBCODE']
control = read_csv('./483OutputA3-3-M-Benchmark.csv', names=NAMES_CSV)
codes_test = set(control.PDBCODE)


def get_control_dataset(code):
    df_control = control[control.PDBCODE == code]
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
    

@pytest.mark.parametrize("code", codes_test)
def test_eval(code):
    df_control = get_control_dataset(code=code)
    df_test = get_test_dataset(code=code)
    testing.assert_frame_equal(df_control, df_test)


# run manually
if __name__ == '__main__':
    code = '6dbp'
    print('2016: \n', get_control_dataset(code))
    print('Today: \n', get_test_dataset(code))