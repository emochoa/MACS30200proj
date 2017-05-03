import math
import numpy as np
import pandas as pd


def value_parser(valuestring):
    '''
    '''

    if valuestring == '-8888' or valuestring == '-9999':
        return np.NaN
    else:
        value = valuestring[0:-1]
        return int(value)


def create_dicto():
    '''
    '''

    names = [str(x) for x in range(1,32)]

    dicto = {}

    for name in names:
        dicto[name] = value_parser

    return dicto




def read_file():
    '''
    '''

    fnm = '/Users/erin/Desktop/persp-data/noaa/ytd-snow-normal.txt'

    n = ['StationID','Month'] + [str(x) for x in range(1,32)]

    d = create_dicto()

    i = ['StationID','Month']

    w = [12,3,10,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7]

    df = pd.read_fwf(fnm, converters = d, names = n, widths = w, index_col = i)

    df = df.mean(axis=1)

    df = df.rename('MeanSnow')

    return df


if __name__ == "__main__":

    df = read_file()

    df.to_csv('snow.csv', header=True)
