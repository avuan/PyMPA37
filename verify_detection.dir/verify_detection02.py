#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# import of useful subroutines and libraries 
# Plot template versus continuous data and check differences

import glob
import json
import math as M
from math import log10
from pprint import pprint

import bottleneck as bn
import matplotlib.pyplot as plt
import numpy as np
import obspy.signal
from obspy.core import *
from obspy.core.event import *
from obspy.core.inventory import read_inventory
from obspy.geodetics import gps2dist_azimuth
from obspy.imaging import *
from obspy.signal import *
from obspy.signal.util import *
from obspy.taup import *
from obspy.taup.taup_create import *


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


def calc_timeshift(eve_lat, eve_lon, eve_dep, sta_lat, sta_lon):
    half_time = 2.5
    epi_dist, az, baz = gps2dist_azimuth(eve_lat, eve_lon, sta_lat, sta_lon)
    epi_dist = epi_dist / 1000
    deg = kilometer2degrees(epi_dist)
    # build_taup_model(taup_model)
    model = TauPyModel(model=taup_model)
    arrivals = model.get_travel_times(source_depth_in_km=eve_dep, distance_in_degree=deg, phase_list=["s", "S"])
    print(arrivals)
    arrS = arrivals[0]
    # correction to be used to evaluate the origin time of event
    origin_time_shift = arrS.time - half_time
    return origin_time_shift


def mag_detect(magt, amaxt, amaxd):
    """
    mag_detect(mag_temp,amax_temp,amax_detect)
    Returns the magnitude of the new detection by using the template/detection amplitude trace ratio
    and the magnitude of the template event
    """
    print("magt == ", magt)
    # print "amaxt == ", amaxt
    # print "amaxd == ", amaxd
    amaxr = amaxt / amaxd
    # print "amaxr == ", amaxr
    magd = magt - log10(amaxr)
    print("magd == ", magd)
    return magd


def read_input_par(vfile):
    with open(vfile) as tf:
        data = tf.read().splitlines()
    stations = data[16].split(" ")
    channels = data[17].split(" ")
    networks = data[18].split(" ")
    lowpassf = float(data[19])
    highpassf = float(data[20])
    tlen_bef = float(data[21])
    tlen_aft = float(data[22])
    UTC_prec = int(data[23])
    cont_dir = "./" + data[24] + "/"
    temp_dir = "./" + data[25] + "/"
    ttimes_dir = "./" + data[26] + "/"
    ev_catalog = str(data[27])
    start_det = int(data[28])
    stop_det = int(data[29])
    det_dur = float(data[30])
    taup_model = str(data[31])
    return stations, channels, networks, lowpassf, highpassf, tlen_bef, tlen_aft, UTC_prec, cont_dir, temp_dir, ttimes_dir, ev_catalog, start_det, stop_det, det_dur, taup_model


def read_sta_inv(invfiles, sta):
    for invv in invfiles:
        inv = read_inventory(invv)
        nt0 = inv[0].select(station=sta)
        if bool(nt0):
            lat = nt0[0].latitude
            lon = nt0[0].longitude
            elev = nt0[0].elevation
            return lat, lon, elev


vfile = 'verify.par'
stations, channels, networks, lowpassf, highpassf, tlen_bef, tlen_aft, UTC_prec, cont_dir, temp_dir, ttimes_dir, ev_catalog, start_det, stop_det, det_dur, taup_model = read_input_par(
    vfile)
invfiles = ["inv.ingv.iv", "inv.ingv.mn"]
Flag_Save_Figure = 1
# pay attention it must be equal to detection

# read Template Catalog
cat = Catalog()
cat = read_events(ev_catalog, format="ZMAP")
ncat = len(cat)
aa = np.loadtxt(ev_catalog)
aa1 = aa[:, 9]
aa2 = aa1 - np.floor(aa1)
aa3 = aa2 * 1000000
print("aa1, aa2, aa3 == ", aa1, aa2, aa3)

# read Catalog as output from Phase Match Filtering
dc = np.loadtxt("./outcat")
# catalog format
# 2013 9 5 20 4 58.74 1.01 0.506 11.818 32 40.3127 0.6997 5.0 5.0
# 2013 9 5 21 43 47.35 1.29 0.394 9.816 105 40.3744 0.702 2.0 4.0
# 2013 9 5 22 48 29.78 0.53 0.516 12.095 111 40.4103 0.6947 2.0 4.0
# 2013 9 6 4 11 10.6 0.35 0.458 11.586 109 40.3846 0.7071 3.0 4.0
# 2013 9 6 5 20 10.96 0.57 0.432 9.892 75 40.4049 0.6818 1.0 5.0

dcyy = dc[:, 0]
dcmm = dc[:, 1]
dcdd = dc[:, 2]
dchh = dc[:, 3]
dcmin = dc[:, 4]
dcss = dc[:, 5]
dcmag = dc[:, 6]
dcavcc = dc[:, 7]
dcthre = dc[:, 8]
dctmp = dc[:, 9]
dclat = dc[:, 10]
dclon = dc[:, 11]
dcdep = dc[:, 12]
dcchan = dc[:, 13]

ndet = len(dc[:, 0])

channel = ["???"]
detection_num = 0

for jf, detection_num in enumerate(range(start_det, stop_det)):
    # for jf, detection_num in enumerate(range(85,86)):
    st_temp = Stream()
    st_cont = Stream()
    st1_temp = Stream()
    st1_cont = Stream()

    # read data from the detection catalog
    yy = int(dcyy[detection_num])
    mm = int(dcmm[detection_num])
    dd = int(dcdd[detection_num])
    hh = int(dchh[detection_num])
    mmn = int(dcmin[detection_num])
    ss = dcss[detection_num]
    print("ss ==", ss)

    detection_otime = UTCDateTime(yy, mm, dd, hh, mmn, ss)
    detection_daily_otime = UTCDateTime(yy, mm, dd, hh, mmn, ss).timestamp - UTCDateTime(yy, mm, dd, 0, 0, 0,
                                                                                         0).timestamp
    magd = dcmag[detection_num]
    avcc = dcavcc[detection_num]
    thre = dcthre[detection_num]
    template_num = int(dctmp[detection_num])
    chan_num = int(dcchan[detection_num])

    # read template data from the template catalog
    eve_lat = cat[template_num].origins[0].latitude
    eve_lon = cat[template_num].origins[0].longitude
    eve_dep = cat[template_num].origins[0].depth / 1000

    # read time and parameters of template event 
    ot = cat[template_num].origins[0].time.datetime
    ot1 = UTCDateTime(ot)
    yyt = ot1.year
    mmt = ot1.month
    ddt = ot1.day
    hht = ot1.hour
    mint = ot1.minute
    sst = ot1.second
    msect = aa3[template_num]
    microsect = int(msect)
    print("aa3[template_num], microsect", aa3[template_num], microsect)

    magt = cat[template_num].magnitudes[0].mag
    template_otime = UTCDateTime(yyt, mmt, ddt, hht, mint, sst, microsect)
    template_daily_otime = UTCDateTime(yyt, mmt, ddt, hht, mint, sst, microsect).timestamp - UTCDateTime(yyt, mmt, ddt,
                                                                                                         0, 0, 0,
                                                                                                         0).timestamp
    print("yyt, mmt, ddt, hht, mint, sst, microsect == ", yyt, mmt, ddt, hht, mint, sst, microsect)
    print("detection_daily_otime, template_daily_otime == ", detection_daily_otime, template_daily_otime)
    # define useful parameters for the stream of continuous data
    # string name
    sday = str(yy)[2:4] + str(mm).zfill(2) + str(dd).zfill(2)
    # filtering range
    bandpass = [lowpassf, highpassf]
    # plotted duration for stream +-det_dur

    # carica template
    stemp_num = str(template_num)
    print("template_num == ", stemp_num)
    for file in glob.glob(temp_dir + stemp_num + '.*'):
        print('file :' + file)
        st_temp += read(file)
    print("Stream ...", st_temp)
    npanels = len(st_temp)
    tt = Trace()
    # calcola il tempo di riferimento minimo dello stream dei template
    print(st_temp)
    refT = min([tt.stats.starttime for tt in st_temp])
    print("refT == ", refT)

    # carica tracce 24h continue
    for file1 in glob.glob(cont_dir + sday + '.*' + '.' + channel[0]):
        print(" file1 ===", file1)
        st_cont += read(file1)
    st_cont.filter('bandpass', freqmin=bandpass[0], freqmax=bandpass[1], zerophase=True)
    # load ttimes calculated by calcTT02.py
    d = {}

    # define travel time file for each template (travel time files for synchronizing CFTs are obtained
    # running calcTT01.py
    travel_file = "%s%s.ttimes" % (ttimes_dir, template_num)

    # store ttimes info in a dictionary
    with open(travel_file, "r") as ttim:
        d = dict(x.rstrip().split(None, 1) for x in ttim)
    ttmin = [float(str(x)) for x in d.values()]
    # print ttmin
    min_time = min(ttmin[:])
    # print d
    # print " min_time == ", min_time

    # select filtered template
    st1_temp = st_temp.select(channel=channel[0])
    tt = Trace()
    tc = Trace()
    count = 0
    nmin = 0

    nfile = len(st1_temp)
    Tstart = np.empty(nfile)
    Tend = np.empty(nfile)
    tdif = np.empty(nfile)
    jf, axarray = plt.subplots(npanels, sharex=True)

    count = 0
    for it, tt in enumerate(st1_temp):
        count = count + 1
        sstat = tt.stats.station
        print("station == ", sstat)
        slat, slon, selev = read_sta_inv(invfiles, sstat)
        print("sta_lat, sta_lon === ", eve_lat, eve_lon, eve_dep, slat, slon)
        ori = calc_timeshift(eve_lat, eve_lon, eve_dep, slat, slon)
        print("ori == ", ori)

        chan = tt.stats.channel
        netwk = tt.stats.network
        idchan = netwk + "." + sstat + ".." + chan
        idchan_dic = netwk + "." + sstat + "." + chan
        print("idchan_dic == ", idchan_dic)
        st1_cont = Stream()
        if st_cont.select(id=idchan).__nonzero__():
            st1_cont = st_cont.select(id=idchan)
            st1_cont.merge()
            print(st1_cont.__str__(extended=True))
            tc = st1_cont[0]
            print("Trace Id.i= ", st1_cont[0])
            tc.trim(starttime=UTCDateTime(detection_otime), endtime=UTCDateTime(detection_otime) + 2 * det_dur,
                    pad=True, nearest_sample=True, fill_value=0)
            print("tc.stats.starttime == ", tc.stats.starttime)
            print("detection_otime == ", UTCDateTime(detection_otime))

            ind_tmin = int((ori) / tc.stats.delta) + 1
            # ind_tmax=int((ori+5.0)/tc.stats.delta)+1
            ind_tmax = int((ori + 10.0) / tc.stats.delta) + 1
            magg = max(abs(tc.data))
            tm = tc.copy()
            amaxat = max(abs(tt.data))
            if tc.data[ind_tmin:ind_tmax].size:
                amaxad = max(abs(tc.data[ind_tmin:ind_tmax]))
            else:
                amaxad = amaxat
                # amaxad=max(abs(tc.data))
            # amaxat=max(abs(tt.data))
            print("amaxat, amaxad == ", amaxat, amaxad)
            amaxmul = amaxad / amaxat
            print("amaxmul == ", amaxmul)
            magd = mag_detect(magt, amaxat, amaxad)
            smagd = str("%4.2f" % magd)
            smagd = "Md=" + smagd + " Mt=" + str(magt)
            mag = amaxad
            tlen = UTCDateTime(detection_otime).timestamp

            # tt.trim(starttime=UTCDateTime(template_otime), endtime=UTCDateTime(ttbegin)+5, fill_value=0)
            ############################################################################################
            ## Plot detected events versus template
            ############################################################################################
            tad = np.arange(0, (tt.stats.npts / tt.stats.sampling_rate), tt.stats.delta)
            t = np.arange(0, tc.stats.npts / tc.stats.sampling_rate, tc.stats.delta)
            axarray[count - 1].plot(t, tc.data, 'k', lw=1.5)
            axarray[count - 1].plot(tad + ori, amaxmul * tt.data, 'r', lw=1.5)
            axarray[count - 1].text(det_dur * 0.05, 0.65 * magg, smagd, fontsize=11)
            axarray[count - 1].text(det_dur * 0.45, 0.65 * magg, 'avcc= ' + str(avcc), fontsize=11)
            axarray[count - 1].text(det_dur * 0.7, 0.65 * magg, 'thre= ' + str(thre), fontsize=11)
            axarray[count - 1].text(det_dur, 0.65 * magg, 'template_num=' + str(template_num), fontsize=11)
            axarray[count - 1].text(det_dur * 1.8, 0.65 * magg, tc.stats.station + '.' + tc.stats.channel, fontsize=11)
            print("det_dur, magg == ", det_dur, magg)
            # axarray[count-1].text(det_dur, -0.9*magg, 'template_time=' + str(tt.stats.starttime), fontsize=11)

            # uncomment 2 lines below to display template and detection time
            # axarray[count-1].text(det_dur, -0.9*magg, 'template_time=' + str(tt.stats.starttime-ori), fontsize=11)
            # axarray[count-1].text(det_dur*0.05, -0.9*magg, 'detection_time=' + str(detection_otime), fontsize=11)
            if count == 1:
                axarray[count - 1].text(det_dur * 0.65, 1.4 * magg,
                                        'Detection = ' + str(detection_num) + ' ' + 'Template = ' + str(template_num),
                                        fontsize=16)
                print("det_dur, magg == ", det_dur, magg)

                # axes = plt.gca()
                # axes.set_ylim([-amaxat,amaxat])
                # axarray[count-1].ylim(-amaxat,amaxat)
            plt.xlabel('Time [s]')
            # plt.ylim(-amaxat,amaxat)
            # fig = plt.gcf()
            # fig.set_size_inches(11.93, 15.98,11.93)
            ###############################################################################
            # Set Flag_Save_Figure=1 to save png files or Flag_Save_Figure=0 to show results
            ###############################################################################
    if Flag_Save_Figure == 0:
        plt.show()
    if Flag_Save_Figure == 1:
        fig = plt.gcf()
        # fig.set_size_inches(11.93,15.98)
        fig.set_size_inches(15.98, 11.93)
        # print "sfig===", sfig
        outfile = sday + sstat + "." + str(detection_num) + "." + str(template_num).zfill(3) + ".png"
        fig.savefig(outfile, dpi=300)
        fig.clf()
