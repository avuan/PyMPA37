#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import numpy as np
import glob
import pandas as pd
import datetime
from obspy import read
from obspy import read_inventory, read_events, Stream, Trace
from obspy.core.utcdatetime import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from obspy.taup.taup_create import build_taup_model
from obspy.taup.tau import TauPyModel


def listdays(year, month, day, period):
    # create a list of days for scanning by templates
    datelist = pd.date_range(datetime.datetime(year, month, day), periods=period).tolist()
    a = list(map(pd.Timestamp.to_pydatetime, datelist))
    days_from_par = []
    for i in a:
        days_from_par.append(i.strftime("%y%m%d"))
    return days_from_par


def create_day_list(catalog, days_from_par):
    # creates day list from the intersection between the datepriod input in trim.par and
    # the template catalog
    days = []
    days_from_catalog = []
    num_eve = len(catalog)
    dd = [catalog[i].origins[0].time.date for i in range(0, num_eve)]
    dd = list(set(dd))

    for day in dd:
        year = str(day.strftime("%y"))
        month = str(day.strftime("%m"))
        dday = str(day.strftime("%d"))
        days_from_catalog += [year + month + dday]

    print(days_from_catalog, days_from_par)

    days = list(sorted(set(days_from_par).intersection(days_from_catalog), key=lambda x:days_from_par.index(x)))

    print(days)
    return days


def kilometer2degrees(kilometer, radius=6371):
    """
    Convenience function to convert kilometers to degrees assuming a perfectly
    spherical Earth.

    :type kilometer: float
    :param kilometer: Distance in kilometers
    :type radius: int, optional
    :param radius: Radius of the Earth used for the calculation.
    :rtype: float
    :return: Distance in degrees as a floating point number.

    """
    pi = 3.14159265359
    return kilometer / (2.0 * radius * pi / 360.0)


def read_input_par(trimfile):

    with open(trimfile) as tf:
        data = tf.read().splitlines()

    stations = data[15].split(" ")
    channels = data[16].split(" ")
    networks = data[17].split(" ")
    lowpassf = float(data[18])
    highpassf = float(data[19])
    tlen_bef = float(data[20])
    tlen_aft = float(data[21])
    utc_prec = int(data[22])
    cont_dir = "./" + data[23] + "/"
    temp_dir = "./" + data[24] + "/"
    dateperiod = data[25].split(" ")
    ev_catalog = str(data[26])
    start_itemp = int(data[27])
    stop_itemp = int(data[28])
    taup_model = str(data[29])

    return (
        stations,
        channels,
        networks,
        lowpassf,
        highpassf,
        tlen_bef,
        tlen_aft,
        utc_prec,
        cont_dir,
        temp_dir,
        dateperiod,
        ev_catalog,
        start_itemp,
        stop_itemp,
        taup_model,
    )


def read_sta_inv(invfile, sta):
    inv = read_inventory(invfile)
    nt0 = inv[0].select(station=sta)
    if nt0:
        lat = nt0[0].latitude
        lon = nt0[0].longitude
        elev = nt0[0].elevation
        print(sta, lat, lon, elev)
    else:
        lat = 999
        lon = 999
        elev = 999
    return lat, lon, elev


trimfile = "./trim.par"


# lat, lon, elev = read_sta_inv(invfile, station)
# print(lat, lon, elev)

[
    stations,
    channels,
    networks,
    lowpassf,
    highpassf,
    tlen_bef,
    tlen_aft,
    utc_prec,
    cont_dir,
    temp_dir,
    dateperiod,
    ev_catalog,
    start_itemp,
    stop_itemp,
    taup_model,
] = read_input_par(trimfile)

# -------
# Define our bandpass min and max values
bandpass = [lowpassf, highpassf]

# Duration to use for template before and after s-waves arrival time in seconds
tmplt_dur = tlen_bef

# ---Read the Catalog of template in zmap_format, filtered by day---#
cat = read_events(ev_catalog, format="ZMAP")
ncat = len(cat)

# ---- The following lines are needed because the input zmap has no decimal year
# ---- in the corresponding column, but fractions of seconds are in the seconds field
aa = np.loadtxt(ev_catalog)
aa1 = aa[:, 9]
aa2 = aa1 - np.floor(aa1)
aa3 = aa2 * 1000000

st = Stream()
st1 = Stream()
st2 = Stream()

# create taup model from a tvel model
tvel = "./" + taup_model + ".tvel"
build_taup_model(tvel)

# generate list of days to process
year = int(dateperiod[0])
month = int(dateperiod[1])
day = int(dateperiod[2])
period = int(dateperiod[3])
days_from_par = listdays(year, month, day, period)
days = create_day_list(cat, days_from_par)

for ista in stations:

    for day in days:
        print("day == ", day)
        inpfiles = cont_dir + day + "." + ista + ".???"
        st.clear()

        for file in glob.glob(inpfiles):
            st += read(file)

        st.merge(method=1, fill_value=0)
        st.detrend("constant")
        st.filter("bandpass", freqmin=bandpass[0], freqmax=bandpass[1], zerophase=True)
        dataYY = int("20" + day[0:2])
        dataMM = int(day[2:4])
        dataDD = int(day[4:6])

        # to avoid errors in the input trim.par at stop_itemp
        if stop_itemp > ncat:
            stop_itemp = ncat

        for iev in range(start_itemp, stop_itemp):
            ot = cat[iev].origins[0].time.datetime
            ot1 = UTCDateTime(ot)
            yy = ot1.year
            mm = ot1.month
            dd = ot1.day
            hh = ot1.hour
            minu = ot1.minute
            sec = ot1.second
            microsec = aa3[iev]
            m = cat[iev].magnitudes[0].mag
            lon = cat[iev].origins[0].longitude
            lat = cat[iev].origins[0].latitude
            dep = cat[iev].origins[0].depth
            # depth in km
            dep = dep / 1000
            microseci = int(microsec)
            ot0 = UTCDateTime(yy, mm, dd, hh, minu, sec, microseci)
            ot2 = UTCDateTime(yy, mm, dd)
            ot3 = UTCDateTime(dataYY, dataMM, dataDD)

            if ot2 == ot3:
                eve_coord = [lat, lon, dep]

            print("ista", ista)

            for network in networks:
                invfile = "./inv." + network
                slat, slon, selev = read_sta_inv(invfile, ista)
                print(ista, slat, slon, selev)

                if slat != 999:
                    print("Station ", ista, " found in inventory ", invfile)
                    break
                else:
                    print("Warning no data found in ", invfile, " for station", ista)
                    continue

            eve_lat = eve_coord[0]
            eve_lon = eve_coord[1]
            eve_dep = eve_coord[2]

            if eve_dep < 1.5:
                eve_dep = 1.5

            epi_dist, az, baz = gps2dist_azimuth(eve_lat, eve_lon, slat, slon)
            epi_dist = epi_dist / 1000
            print("epi_dist==", epi_dist)
            deg = kilometer2degrees(epi_dist)
            print("deg==", deg)
            print("eve_dep==", eve_dep)
            model = TauPyModel(model=taup_model)
            arrivals = model.get_travel_times(
                source_depth_in_km=eve_dep,
                distance_in_degree=deg,
                phase_list=["s", "S"],
            )
            arrS = arrivals[0]
            print("arrS.time=...", arrS.time)

            stime = UTCDateTime(ot0) + arrS.time - tlen_bef
            print("stime", stime)
            etime = UTCDateTime(ot0) + arrS.time + tlen_aft
            print("etime", etime)

            # cut the 3-component template and save file
            nchannels = len(channels)

            for ichan in range(0, nchannels):
                print("ista", ista)
                st1.clear()
                # print("FILE", file)
                st1 = st.copy()
                tw = Trace()
                st2.clear()
                print(st1.select(station=ista, channel=channels[ichan]))
                st2 = st1.select(station=ista, channel=channels[ichan])

                if st2.__nonzero__():
                    tw = st2[0]

                    if tw.trim(stime, etime).__nonzero__():
                        print(tw)
                        netwk = tw.stats.network
                        ch = tw.stats.channel
                        tw.trim(stime, etime)
                        newfile = (
                            temp_dir
                            + str(iev)
                            + "."
                            + netwk
                            + "."
                            + ista
                            + ".."
                            + ch
                            + ".mseed"
                        )
                        print(newfile)
                        tw.write(newfile, format="MSEED")
                    else:
                        pass

                else:
                    pass
