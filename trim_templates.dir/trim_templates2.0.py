#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import glob
import json
import obspy.signal
import numpy as np
import matplotlib.pyplot as plt

from math import log10
from obspy.core import *
from obspy.signal import *
from obspy.signal.util import *
from obspy.imaging import *
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import *
from obspy.taup.taup_create import *
from obspy.io.mseed.util import get_timing_and_data_quality
from obspy import Catalog, UTCDateTime
from obspy.core.event import * 
from obspy.io.zmap import *

## import end -----------------------------------------------------------

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
    pi=3.14159265359
    return kilometer / (2.0 * radius * pi / 360.0)

#-------
# Define our bandpass min and max values
bandpass=[2.0,8.0] 

# Percent taper to apply to the template
FlagTaper=0
taper=0.1

# Duration to use for template before and after s-waves arrival time in seconds       
tmplt_dur=2.5

#---Read the Catalog of template in zmap_format, filtered by day---#

cat = Catalog()
cat = read_events("templates.zmap", format="ZMAP")
ncat = len(cat)

aa = np.loadtxt("templates.zmap")
aa1 = aa[:,9]
aa2 = aa1 - np.floor(aa1)
aa3 = aa2 * 1000000
#print(cat.__str__(print_all=True))

# setup station list
stations=["AQU", "CAMP", "CERT", "FAGN", "FIAM", "GUAR", "INTR", "MNS", "NRCA", "TERO"]

# setup station dictionary, no need to order it
sta_coord={'AQU': [42.35400, 13.40500, 0.7], 'CAMP': [42.53578, 13.40900, 1.3], \
    'CERT': [41.94903, 12.98176, 0.773], 'FAGN': [42.26573, 13.58379, 0.8], \
    'FIAM': [42.26802, 13.11718, 1.070], 'GUAR': [41.79450, 13.31229, 0.7], \
    'INTR': [42.01154, 13.90460, 0.924], 'MNS': [42.38546, 12.68106, 0.706], \
    'NRCA': [42.83355, 13.11427, 0.927], 'TERO': [42.62279, 13.60393, 0.673]}


st = Stream()
st1 = Stream()
st2 = Stream()
flist = "./lista1"
fname="%s" % (flist)

# array of days is built deleting last line character (/newline) ls -1 command include a newline character at the end
with open(fname) as fl:
    days = [line[:-1] for line in fl]
    print(days)
fl.close()

for day in days:
    print("day == ", day)
    inpfiles="./24h/" + day + ".*.???"
    st.clear()
    for file in glob.glob(inpfiles):
        st+=read(file)
    st.merge(method = 1, fill_value = 0)
    st.detrend('constant')
    st.filter('bandpass', freqmin = bandpass[0], freqmax = bandpass[1], zerophase = True)
    dataYY=int("20" + day[0:2])
    dataMM=int(day[2:4])
    dataDD=int(day[4:6])
    #for iev in range(0,ncat):
    for iev in range(26, 27):
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
        dep = dep/1000
        microseci = int(microsec)
        ot0 = UTCDateTime(yy, mm, dd, hh, minu, sec, microseci)
        ot2 = UTCDateTime(yy, mm, dd)
        ot3 = UTCDateTime(dataYY, dataMM, dataDD)
      
        if ot2==ot3:
            eve_coord=[lat, lon, dep]
            channels=["BHE", "BHN", "BHZ"]
            for ista in stations:
                print("ista", ista)
                sta_lat = sta_coord.get(ista)[0]
                sta_lon = sta_coord.get(ista)[1]
                eve_lat = eve_coord[0]
                eve_lon = eve_coord[1]
                eve_dep = eve_coord[2]
                if eve_dep < 1.5:
                    eve_dep = 1.5
                epi_dist, az, baz = gps2dist_azimuth(eve_lat, eve_lon, sta_lat, sta_lon)
                epi_dist = epi_dist / 1000
                print("epi_dist==", epi_dist)
                deg=kilometer2degrees(epi_dist)
                print("deg==", deg)
                print("eve_dep==", eve_dep)
                model = TauPyModel(model="aquila_kato")
                arrivals = model.get_travel_times(source_depth_in_km=eve_dep, distance_in_degree=deg, phase_list=["s", "S"])
                arrS = arrivals[0]
                print("arrS.time=...", arrS.time)
							
                stime=UTCDateTime(ot0)+arrS.time-tmplt_dur
                print("stime", stime) 
                etime=UTCDateTime(ot0)+arrS.time+3*tmplt_dur
                print("etime", etime)
	
                # cut the 3-component template and save file
                nchannels=len(channels)
                for ichan in range(0,nchannels):
                    print("ista", ista)
                    st1.clear()
                    #print("FILE", file)
                    st1=st.copy()
                    #tq = get_timing_and_data_quality(file)
                    #for k, v in sorted(tq.items()):
                    #     print(k, v)
                    tw=Trace()
                    st2.clear()
                    print(st1.select(station=ista, channel=channels[ichan]))
                    st2 = st1.select(station=ista, channel=channels[ichan])
                    if st2.__nonzero__():
                        tw=st2[0]
                        #for k, v in sorted(st2[0].stats.mseed.items()):
                        #    print("'%s': %s" % (k, str(v))) 
#                        print("delta == ", tw.stats.delta)
#                       print("starttime ==", tw.stats.starttime)
#                        print("endtime == ", tw.stats.endtime)
#                        print("calib ==", tw.stats.calib)
#                        print("sampling_rate == ", tw.stats.sampling_rate)
#                        print("network == ", tw.stats.network)
#                        print("channel == ", tw.stats.channel)
#                        print("station == ", tw.stats.station)
#                        print("tw.data == ", tw.data)
                        if tw.trim(stime,etime).__nonzero__():
                            print(tw)
                            netwk = tw.stats.network
                            ch=tw.stats.channel
                            tw.trim(stime,etime)
	
                            #tw.plot()
                            newfile="./template/" + str(iev) + "." + netwk + "." + ista + ".." + ch + ".mseed"
                            #print ot0
                            print(newfile)
                            tw.write(newfile, format = "MSEED")
                        else:
                            pass
                    else:
                        pass
