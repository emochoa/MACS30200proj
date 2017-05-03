import math
import numpy as np
import pandas as pd
from datetime import date as date, datetime as dt
from dateutil.relativedelta import relativedelta as rd

def date_fixer(datestring):
    '''
    Fixes dates that come in 'YYYY-MM-DD' format and changes them to datetime.
    Inputs:
        date_int (string): A date in 'YYYY-MM-DD' format.
    Output:
        None.
    Returns:
        date (datetime object): A corrected date as a datetime object (or
                                np.NaN if the date was invalid).
    '''

    if datestring:
        date = dt.strptime(datestring,'%Y-%m-%d')
        return date
    else:
        return np.NaN


def read_file():
    '''
    '''

    fnm = '/Users/erin/Desktop/persp-data/sunshine.csv'

    n = ['StationID','Date','DayLengthMins','SolarNoonUTC']

    c = ['StationID','Date','DayLengthMins']

    d = {'Date': date_fixer}

    df = pd.read_csv(fnm, usecols = c, converters = d, names = n, skiprows = 1)

    df['Year'] = df.Date.apply(lambda x: x.year)
    df['Month'] = df.Date.apply(lambda x: x.month)

    df = df.groupby(['StationID','Year','Month']).mean()

    return df


if __name__ == "__main__":

    df = read_file()

    df.to_csv('mean_daylight_minutes.csv')
