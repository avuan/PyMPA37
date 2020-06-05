#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

# To run calcTT05.py template.zmat (template catalog) is needed
# Synchronization of CFTs is based on delays in cutting templates
# Origin time is retrieved from travel time of the epicenter closest station
# minus the half_length in seconds of template events
# The program runs before computing CFTs and stacking to obtain CCMAD and STD

import glob
import os
import numpy as np
from obspy.core import read, Stream, Trace
from obspy.core.event import read_events
from obspy.core.inventory import read_inventory
from obspy.geodetics import gps2dist_azimuth
from obspy.taup.tau import TauPyModel
from obspy.io.quakeml.core import _is_quakeml
from obspy.io.zmap.core import _is_zmap


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


def read_input_par(filepar):
    with open(filepar) as tf:
        data = tf.read().splitlines()
    stations = data[14].split(" ")
    channels = data[15].split(" ")
    networks = data[16].split(" ")
    lowpassf = float(data[17])
    highpassf = float(data[18])
    tlen_bef = float(data[19])
    tlen_aft = float(data[20])
    utc_prec = int(data[21])
    tempinp_dir = "./" + data[22] + "/"
    tempout_dir = "./" + data[23] + "/"
    ev_catalog = str(data[24])
    start_itemp = int(data[25])
    stop_itemp = int(data[26])
    taup_model = str(data[27])
    return (
        stations,
        channels,
        networks,
        lowpassf,
        highpassf,
        tlen_bef,
        tlen_aft,
        utc_prec,
        tempinp_dir,
        tempout_dir,
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


filepar = "./times.par"

[
    stations,
    channels,
    networks,
    lowpassf,
    highpassf,
    tlen_bef,
    tlen_aft,
    utc_prec,
    tempinp_dir,
    tempout_dir,
    ev_catalog,
    start_itemp,
    stop_itemp,
    taup_model,
] = read_input_par(filepar)
st = Stream()
tr = Trace()

# half time length of templates is used to calculate origin_time
half_time = tlen_bef

# read event coordinates from catalog
print("event catalog should be ZMAP or QUAKEML")

if _is_zmap(ev_catalog):
    print("reading ZMAP catalog")

elif _is_quakeml(ev_catalog):
    print("reading QUAKEML catalog")

else:
    print("warning error in reading ZMAP or QUAKEML")

cat = read_events(ev_catalog)

ncat = len(cat)
print(cat.__str__(print_all=True))

for iev in range(start_itemp, stop_itemp):
    inplist = tempinp_dir + str(iev) + ".??" + "." + "*" + "..???" + "." + "mseed"
    print("inplist == ...", inplist)
    st.clear()

    for filename in glob.glob(inplist):
        csize = os.path.getsize(filename)

        if csize > 0:
            print("filename ..", filename)
            st += read(filename)

    nst = len(st)
    # print nst
    print(st.__str__(extended=True))
    refT = min([tr.stats.starttime for tr in st])
    Tshift = np.empty(nst)
    arrP = np.empty(nst)
    arrS = np.empty(nst)
    origin_time_shift = np.empty(1)
    fout1 = tempout_dir + str(iev) + ".ttimes"
    fileout = open(fout1, "w+")

    for it, tr in enumerate(st):
        Tshift[it] = tr.stats.starttime - refT

        if Tshift[it] == 0:
            # define epicenter closest station
            ref_sta = tr.stats.station
            # event details
            ot = cat[iev].origins[0].time
            m = cat[iev].magnitudes[0].mag
            lon = cat[iev].origins[0].longitude
            lat = cat[iev].origins[0].latitude
            dep = cat[iev].origins[0].depth
            # depth in km
            dep = dep / 1000
            print("time, mag, lon, lat, dep", ot, m, lon, lat, dep)
            eve_coord = [lat, lon, dep]

            for network in networks:
                invfile = "./inv." + network
                slat, slon, selev = read_sta_inv(invfile, ref_sta)
                print(ref_sta, slat, slon, selev)

                if slat != 999:
                    print("Station ", ref_sta, " found in inventory ", invfile)
                    break
                else:
                    print("Warning no data found in ", invfile, " for station", ref_sta)
                    continue

            eve_lat = eve_coord[0]
            eve_lon = eve_coord[1]
            eve_dep = eve_coord[2]
            # print "sta_lon, sta_lat, eve_lon, eve_lat==", sta_lon,
            # sta_lat, eve_lon, eve_lat
            epi_dist, az, baz = gps2dist_azimuth(eve_lat, eve_lon, slat, slon)
            epi_dist = epi_dist / 1000
            deg = kilometer2degrees(epi_dist)
            model = TauPyModel(model=taup_model)
            arrivals = model.get_travel_times(
                source_depth_in_km=eve_dep,
                distance_in_degree=deg,
                phase_list=["s", "S"],
            )
            print(arrivals)
            # arrP = arrivals[0]
            arrS = arrivals[0]
            # correction to be used to evaluate the origin time of event
            origin_time_shift = arrS.time - half_time
            print("OT...= ", origin_time_shift)

    for ii, tr in enumerate(st):
        Tshift[ii] = Tshift[ii] - origin_time_shift
        ssta = tr.stats.station
        schan = tr.stats.channel
        snet = tr.stats.network
        s1 = (
            str(snet)
            + "."
            + str(ssta)
            + "."
            + str(schan)
            + " "
            + str(Tshift[ii])
            + "\n"
        )
        fileout.write(s1)

    fileout.close()
