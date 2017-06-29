#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
###AV2015
# To run calcTT01.py template.zmat (template catalog) is needed
# Synchronization of CFTs is based on delays in cutting templates
# Origin time is retrieved from travel time of the epicenter closest station
# minus the half_length in seconds of template events
# The program runs before computing CFTs and stacking to obtain CCMAD and STD

import glob
import json

from math import log10
import obspy.signal
from obspy.core import *
from obspy.signal import *
from obspy.signal.util import *
from obspy.io.xseed import Parser
from obspy.io.xseed import DEFAULT_XSEED_VERSION, utils, blockette
from obspy.io.xseed.utils import SEEDParserException
from obspy.imaging import *
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import *
from obspy.taup.taup_create import *
from obspy import Catalog, UTCDateTime
from obspy.core.event import *
from obspy.io.zmap import *
from mpl_toolkits.basemap import Basemap
import numpy as np
import bottleneck as bn
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from pprint import pprint


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

    .. rubric:: Example

    >>> from obspy.core.util import kilometer2degrees
    >>> kilometer2degrees(300)
    2.6979648177561915
    """
    pi = 3.14159265359
    return kilometer / (2.0 * radius * pi / 360.0)


st = Stream()
tr = Trace()
##### important hal_time variable #######
# half time length of templates is used to calculate origin_time
half_time = 2.5

# defines directories of template events)
temp_dir = "./template/"
foutdir = "./ttimes/"

# read event coordinates from catalog
cat = Catalog()
cat = read_events("templates.zmap", format="ZMAP")
# cat.write("my_cat.xml", format="QUAKEML")
# cat.plot(projection="local", resolution= "l", color="date")
ncat = len(cat)
# print(cat[0])
# print(cat[1])
# print(cat[511])
print(cat.__str__(print_all=True))

netwk = ["MN", "IV"]
# setup station list
stations = ["AQU", "CAMP", "CERT", "FAGN", "FIAM", "GUAR", "INTR", "MNS", "NRCA", "TERO"]

# setup station dictionary, no need to order it
d1 = {'AQU': [42.35400, 13.40500, 0.7], 'CAMP': [42.53578, 13.40900, 1.3], \
      'CERT': [41.94903, 12.98176, 0.773], 'FAGN': [42.26573, 13.58379, 0.8], \
      'FIAM': [42.26802, 13.11718, 1.070], 'GUAR': [41.79450, 13.31229, 0.7], \
      'INTR': [42.01154, 13.90460, 0.924], 'MNS': [42.38546, 12.68106, 0.706], \
      'NRCA': [42.83355, 13.11427, 0.927], 'TERO': [42.62279, 13.60393, 0.673]}

for iev in range(0, ncat + 1):
    # for iev in range(25,26):

    inplist = temp_dir + str(iev) + ".??" + "." + "*" + "..???" + "." + "mseed"
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
    fout1 = foutdir + str(iev) + ".ttimes"
    fileout = open(fout1, 'w+')
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
            #
            sta_lat = d1.get(ref_sta)[0]
            sta_lon = d1.get(ref_sta)[1]
            eve_lat = eve_coord[0]
            eve_lon = eve_coord[1]
            eve_dep = eve_coord[2]
            # print "sta_lon, sta_lat, eve_lon, eve_lat==", sta_lon, sta_lat, eve_lon, eve_lat
            epi_dist, az, baz = gps2dist_azimuth(eve_lat, eve_lon, sta_lat, sta_lon)
            epi_dist = epi_dist / 1000
            deg = kilometer2degrees(epi_dist)
            model = TauPyModel(model="aquila_kato")
            arrivals = model.get_travel_times(source_depth_in_km=eve_dep, distance_in_degree=deg, phase_list=["s", "S"])
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
        s1 = str(snet) + "." + str(ssta) + "." + str(schan) + " " + str(Tshift[ii]) + "\n"
        fileout.write(s1)
    fileout.close()
