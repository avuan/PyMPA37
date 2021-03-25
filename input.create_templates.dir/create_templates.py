#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import glob
import datetime
from obspy import read
from time import perf_counter
from obspy import read_inventory, read_events, Stream, Trace
from obspy.core.utcdatetime import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from obspy.taup.taup_create import build_taup_model
from obspy.taup.tau import TauPyModel
from obspy.io.zmap.core import _is_zmap
from obspy.io.quakeml.core import _is_quakeml


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

    stations = data[16].split(" ")
    channels = data[17].split(" ")
    networks = data[18].split(" ")
    lowpassf = float(data[19])
    highpassf = float(data[20])
    tlen_bef = float(data[21])
    tlen_aft = float(data[22])
    utc_prec = int(data[23])
    time_length_window = int(data[24])
    cont_dir = "./" + data[25] + "/"
    temp_dir = "./" + data[26] + "/"
    ev_catalog = str(data[27])
    start_itemp = int(data[28])
    stop_itemp = int(data[29])
    taup_model = str(data[30])
    limit_epi_dist = float(data[31])

    return (
        stations,
        channels,
        networks,
        lowpassf,
        highpassf,
        tlen_bef,
        tlen_aft,
        utc_prec,
        time_length_window,
        cont_dir,
        temp_dir,
        ev_catalog,
        start_itemp,
        stop_itemp,
        taup_model,
        limit_epi_dist,
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

start_time = perf_counter()
trimfile = "./trim.par"


[
    stations,
    channels,
    networks,
    lowpassf,
    highpassf,
    tlen_bef,
    tlen_aft,
    utc_prec,
    time_length_window,
    cont_dir,
    temp_dir,
    ev_catalog,
    start_itemp,
    stop_itemp,
    taup_model,
    limit_epi_dist,
] = read_input_par(trimfile)

# -------
# Define our bandpass min and max values
bandpass = [lowpassf, highpassf]

# Duration to use for template before and after s-waves arrival time in seconds
tmplt_dur = tlen_bef

# ---Read the Catalog of template in zmap_format, filtered by day---#
# cat = read_events(ev_catalog, format="ZMAP")
cat = read_events(ev_catalog)
ncat = len(cat)

# ---- The following lines are needed because the input zmap has no decimal year
# ---- in the corresponding column, but fractions of seconds are in the seconds field

# check of catalog file if is zmap take microseconds from the last column

if _is_zmap(ev_catalog):
    print("reading ZMAP catalog")
    aa = np.loadtxt(ev_catalog)

    if ncat > 1:
        aa1 = aa[:, 9]
    elif ncat == 1:
        aa1 = aa[9]
    aa2 = aa1 - np.floor(aa1)
    aa3 = aa2 * 1000000

elif _is_quakeml(ev_catalog):
    print("reading QUAKEML catalog")

else:
    print("warning error in reading ZMAP or QUAKEML")

st = Stream()
st1 = Stream()
st2 = Stream()

# create taup model from a tvel model
tvel = "./" + taup_model + ".tvel"
build_taup_model(tvel)

for iev in range(start_itemp, stop_itemp):
    ot = cat[iev].origins[0].time.datetime
    ot1 = UTCDateTime(ot)
    yy = ot1.year
    mm = ot1.month
    dd = ot1.day
    hh = ot1.hour
    minu = ot1.minute
    sec = ot1.second

    if _is_zmap(ev_catalog):

        if ncat == 1:
            microsec = aa3

        else:
            microsec = aa3[iev]

    else:
        microsec = ot1.microsecond

    microsec = int(microsec)
    m = cat[iev].magnitudes[0].mag
    lon = cat[iev].origins[0].longitude
    lat = cat[iev].origins[0].latitude
    dep = cat[iev].origins[0].depth / 1000
    eve_coord = [lat, lon, dep]
    print(eve_coord)

    ot0 = UTCDateTime(yy, mm, dd, hh, minu, sec, microsec)

    # set startime and stoptime values to reduce memory charge
    ot_start = UTCDateTime(yy, mm, dd, hh, minu, sec, microsec) - time_length_window
    ot_stop = UTCDateTime(yy, mm, dd, hh, minu, sec, microsec) + time_length_window

    day = str(yy)[2:4] + str(mm).zfill(2) + str(dd).zfill(2)
    print("day == ", day)
    inpfiles = cont_dir + day + "." + "*.???"

    st.clear()

    for file in glob.glob(inpfiles):
        st += read(file, starttime=ot_start, endtime=ot_stop)

    st.merge(method=1, fill_value=0)
    st.detrend("constant")
    st.filter("bandpass", freqmin=bandpass[0], freqmax=bandpass[1], zerophase=True)
    print(st)

    for tr in st:

        ista = tr.stats.station
        network = tr.stats.network
        channel = tr.stats.channel
        id = tr.id
        # to avoid errors in the input trim.par at stop_itemp
        print("ista, network == ", ista, network)

        invfile = "./inv." + network
        print(id)
        slat, slon, selev = read_sta_inv(invfile, ista)
        print(ista, slat, slon, selev)

        if slat != 999:
            print("Station ", ista, " found in inventory ", invfile)
            # break
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

        if epi_dist < limit_epi_dist:
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

            if tr.trim(stime, etime).__nonzero__():
                print(tr)
                netwk = tr.stats.network
                ch = tr.stats.channel
                tr.trim(stime, etime)
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
                print(tr.stats.mseed.encoding)
                tr.write(newfile, format="MSEED")
            else:
                pass
        else:
            pass
    else:
        pass
print(" elapsed time ", perf_counter() - start_time, " seconds")
