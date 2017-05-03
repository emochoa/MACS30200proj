import math
import numpy as np
import pandas as pd


GRP = ['StationID','Month']

def value_parser(valuestring):
    '''
    '''

    if valuestring == '-8888' or valuestring == '-9999':
        return np.NaN
    else:
        value = valuestring[0:-1]
        return int(value) / 10.0


def create_dicto():
    '''
    '''

    names = [str(x) for x in range(1,32)]

    dicto = {}

    for name in names:
        dicto[name] = value_parser

    return dicto


def monthly_mean(df,measure):
    '''
    '''
    monthly_mean_df = df.mean(axis=1).reset_index().groupby(GRP).mean()
    monthly_mean_df = monthly_mean_df.rename(index=str,columns={0:measure})

    return monthly_mean_df


def monthly_min(df,measure):
    '''
    '''
    monthly_min_df = df.min(axis=1).reset_index().groupby(GRP).min()
    monthly_min_df = monthly_min_df.rename(index=str, columns={0:measure})

    return monthly_min_df

def monthly_max(df,measure):
    '''
    '''
    monthly_max_df = df.max(axis=1).reset_index().groupby(GRP).max()
    monthly_max_df = monthly_max_df.rename(index=str,columns={0:measure})

    return monthly_max_df

def read_file(dicto):
    '''
    '''

    feature = dicto['feature']
    fnm = dicto['fname']
    measures = dicto['measures']
    fnm_out = dicto['fname_out']

    i = GRP

    n = i + ['Day'] + [str(x) for x in range(1,25)]

    c = n[0:2] + n[3:]

    d = create_dicto()

    w = [11,3,3] + [7] * 24

    df = pd.read_fwf(fnm, converters = d, names = n, widths = w, index_col = i,
                     usecols = c)

    monthly_mean_df = monthly_mean(df,measures[0])

    df_list = []

    #'WindChill','Temperature','Pressure','HeatIndex','DewPoint'
    if feature != 'HeatIndex':
        monthly_min_df = monthly_min(df,measures[-1])
        df_list.append(monthly_min_df)

    if feature != 'WindChill':
        monthly_max_df = monthly_max(df,measures[1])
        df_list.append(monthly_max_df)

    df2 = monthly_mean_df.join(df_list)

    write_file(df2, fnm_out)


def write_file(df, fname):
    '''
    '''

    df.to_csv(fname, header=True)


def create_dictionaries():
    '''
    '''

    path = '/Users/erin/Desktop/persp-data/noaa/'

    features = ['WindChill','Temperature','Pressure','HeatIndex','DewPoint']

    fnames = ['hly-wchl-normal','hly-temp-normal','hly-pres-normal',
              'hly-hidx-normal','hly-dewp-normal']

    measures = [tuple(('MonthlyMeanWindChill','MinMeanWindChill')),
                tuple(('MonthlyMeanTemp','MaxMeanTemp','MinMeanTemp')),
                tuple(('MonthlyMeanPressure','MaxMeanPressure','MinMeanPressure')),
                tuple(('MonthlyMeanHeatIdx','MaxMeanHeadIdx')),
                tuple(('MonthlyMeanDewPt','MaxMeanDewPt','MinMeanDewPt'))]

    dicto_list = []

    for i, feature in enumerate(features):
        dicto = {}
        dicto['feature'] = feature
        dicto['fname'] = path + fnames[i] + '.txt'
        dicto['fname_out'] = path +fnames[i] + '.csv'
        dicto['measures'] = measures[i]
        dicto_list.append(dicto)

    return dicto_list



if __name__ == "__main__":

    dicto_list = create_dictionaries()

    for dicto in dicto_list:
        read_file(dicto)
