import math
import astral
import numpy as np
import pandas as pd
from datetime import date, datetime
from dateutil.relativedelta import relativedelta

FNAME = 'sunshine.csv'
HEADERS = 'StationID,Date,DayLengthMins,SolarNoonUTC\n'

def build_lat_lon_list():
    '''
    '''

    fname = '/Users/erin/Desktop/persp-data/noaa/allstations.txt'

    c = ['ID','LATITUDE','LONGITUDE','ELEVATION','STATE']

    n = ['ID','LATITUDE','LONGITUDE','ELEVATION','STATE','NAME','GSNFLAG',
          'HCNFLAG','WMOID','METHOD']

    w = [11,10,10,7,2,32,3,4,1,5]

    df = pd.read_fwf(fname, header = None, usecols = c, names = n, widths = w)


    #DROP STATES/TERRITORIES NOT IN THE CONTIGUOUS 48
    df = df[df.STATE != 'AK']
    df = df[df.STATE != 'HI']
    df = df[df.STATE != 'AS']
    df = df[df.STATE != 'ON']
    df = df[df.STATE != 'MP']
    df = df[df.STATE != 'FM']
    df = df[df.STATE != 'GU']
    df = df[df.STATE != 'UM']
    df = df[df.STATE != 'PW']
    df = df[df.STATE != 'MH']
    df = df[df.STATE != 'PR']
    df = df[df.STATE != 'VI']

    #df['TimeZone'] = TIMEZONE[df.State]

    ids = list(df.ID)
    lat = list(df.LATITUDE)
    lon = list(df.LONGITUDE)
    #state = list(df.STATE)
    #tz = list(df.TimeZone)
    #elev = list(df.ELEVATION)

    table = []

    for i in range(len(ids)):
        #table.append(tuple((ids[i],lat[i],lon[i],state[i],tz[i],elev[i])))
        table.append(tuple((ids[i],lat[i],lon[i])))

    return table


def build_date_range():
    '''
    '''

    bgn = date(1981,1,1)
    end = date(2010,12,31)

    datelist = []

    dt = bgn

    #Every date in the period
    while dt <= end:
        datelist.append(dt)
        dt += relativedelta(days=1)

    '''
    # First, last, and ides
    while dt <= end:
        next_month = dt + relativedelta(months=1)
        end_of_month = next_month + relativedelta(days=-1)
        ides = dt + relativedelta(days = (math.floor(end_of_month.day / 2) - 1))
        datelist.append(dt)
        datelist.append(ides)
        datelist.append(end_of_month)
        dt += relativedelta(months=1)
    '''

    return datelist


def write_headers():
    '''
    '''

    sunshine = open(FNAME, 'w')

    sunshine.write(HEADERS)

    sunshine.close()


def write_data(data_list):
    '''
    '''

    sunshine = open(FNAME, 'a')

    for line in data_list:
        sunshine.write(line)

    sunshine.close()


def build_daylength(station_record):
    '''
    '''

    l = astral.Location()

    l.name, l.latitude, l.longitude = station_record

    data_list = []
    dates = build_date_range()


    for date in dates:
        try:
            sun = l.sun(date)
            sunrise = sun['sunrise']
            sunset = sun['sunset']
            noon = sun['noon']
            daylength = sunset - sunrise
            daylight_mins = daylength.seconds / 60

            data = l.name + ',' + str(date) + ',' + str(daylight_mins) + \
                   ',' + str(noon) + '\n'

            data_list.append(data)

        except:
            print(l.name, str(date))
            continue

    write_data(data_list)



if __name__ == "__main__":

    write_headers()

    table = build_lat_lon_list()

    for station_record in table:
        build_daylength(station_record)
