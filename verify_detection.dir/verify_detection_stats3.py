#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# import of useful subroutines and libraries 
# Plot template versus continuous data and check differences

import glob
import json
from math import log10
import obspy.signal
from obspy.core import *
from obspy.core.event import *
from obspy.signal import *
from obspy.signal.util import *
from obspy.imaging import *
import numpy as np
import math as M
import bottleneck as bn
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from pprint import pprint
from matplotlib.backends.backend_pdf import PdfPages
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import *
from obspy.taup.taup_create import *
import os.path


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
    # build_taup_model('/Users/msugan/work/italia_giappone_2015/test_py/test7/aquila_kato.tvel')
    model = TauPyModel(model="aquila_kato")
    arrivals = model.get_travel_times(source_depth_in_km=eve_dep, distance_in_degree=deg, phase_list=["s", "S"])
    # print arrivals
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
    # print "magt == ", magt
    amaxr = amaxt / amaxd
    # print "amaxr == ", amaxr
    magd = magt - log10(amaxr)
    # print "magd == ", magd
    return magd


# stations=["ALCN", "ALCX", "CMAS", "COBS"]
# stations=["ALCN"]
# channels=["EHE", "EHN", "EHZ"]

netwk = ["IV", "MN"]
sta_coord = {'ORI': [40.05096, 16.450405, 0.375], 'SALB': [39.8772, 16.3459, 1.200], 'SCHR': [40.19924, 16.0759, 0.968], \
             'SIRI': [40.1821, 15.8675, 1.063], 'CUC': [39.9938, 15.8155, 0.637], 'MMN': [39.890961, 15.990414, 0.921], \
             'T0701': [39.98617, 16.116141, 0.882], 'T0711': [39.93571, 16.06183, 0.751],
             'T0715': [39.8384, 16.0683, 1.044], 'CET2': [39.528756, 15.954618, 0.675]}


UTCDateTime.DEFAULT_PRECISION = 6
Flag_Save_Figure = 1
# pay attention it must be equal to detection

# read Template Catalog
cat = Catalog()
cat = read_events("./templates.zmap", format="ZMAP")
ncat = len(cat)
aa = np.loadtxt("./templates.zmap")
aa1 = aa[:, 9]
aa2 = aa1 - np.floor(aa1)
aa3 = aa2 * 1000000
# print "aa1, aa2, aa3 == ", aa1, aa2, aa3

# read Catalog as output from Phase Match Filtering
# dc =  np.loadtxt("./outcat")
# dc =  np.loadtxt("./emilia4_detected_12_3_min25.txt")
# dc =  np.loadtxt("/Volumes/Amatrice/Amatrice_data/output/figure/outcat.key.selected.sort_12_3_03.txt")
# dc =  np.loadtxt("/Volumes/Amatrice/Amatrice_data/amatrice_trim/outcat_3_3_0.5.no_match-key_verificati.format")
# dc =  np.loadtxt("/Volumes/Amatrice/Amatrice_data/amatrice_trim/key_verificati-outcat_3_3_0.5.no_match.format")
dc = np.loadtxt("./outcat")

# catalog format
# 2013 9 5 20 4 58.74 1.01 0.506 11.818 32 40.3127 0.6997 5.0 5.0


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

# for jf, detection_num in enumerate(range(2450,2451)):
for jf, detection_num in enumerate(range(0, ndet)):
    # for jf, detection_num in enumerate(range(600,800)):
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
    # print "ss ==", ss

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

    # read stats from file *stats
    yymmdd = str(yy)[2:4] + str(mm).zfill(2) + str(dd).zfill(2)
    stats_dir = "./"
    filestats = "%s%s.%s.stats" % (stats_dir, str(int(template_num)), yymmdd)
    stat = open(filestats, 'r')
    times = str(detection_otime)
    linestart = 0
    linestop = 0
    istop = 0
    for ist, lstat in enumerate(stat):
        lstat = lstat.strip()
        # print("lstat == ", lstat)
        # print("times == ", times)
        if yymmdd in lstat and times not in lstat and istop == 0:
            linestart = ist
        if yymmdd in lstat and times in lstat:
            linestop = ist
            istop = 1
        # print("1 linestart, linestop == ", linestart, linestop)
    stat.close()
    stat = open(filestats, 'r')
    lines = []
    for ichan, lstat1 in enumerate(stat):
        lstat1 = lstat1.strip()
        # print("2 linestart, linestop == ", linestart, linestop)
        if ichan >= linestart and ichan < linestop:
            # print("lstat1 == ", lstat1)
            cols = lstat1.split()
            channel_name = cols[0]
            channel_cmp = cols[1]
            channel_ccros = cols[3]
            channel_nsample = cols[4]
            lines.append(cols)

        # indices = [i for i, s in enumerate(lines) if 'MN.AQU'in s and 'HHN' in s]

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
    # print "aa3[template_num], microsect", aa3[template_num], microsect

    magt = cat[template_num].magnitudes[0].mag
    template_otime = UTCDateTime(yyt, mmt, ddt, hht, mint, sst, microsect)
    template_daily_otime = UTCDateTime(yyt, mmt, ddt, hht, mint, sst, microsect).timestamp - UTCDateTime(yyt, mmt, ddt,
                                                                                                         0, 0, 0,
                                                                                                         0).timestamp
    # print "yyt, mmt, ddt, hht, mint, sst, microsect == ", yyt, mmt, ddt, hht, mint, sst, microsect
    # print "detection_daily_otime, template_daily_otime == ", detection_daily_otime, template_daily_otime
    # define useful parameters for the stream of continuous data
    # string name
    sday = str(yy)[2:4] + str(mm).zfill(2) + str(dd).zfill(2)
    # filtering range
    # bandpass=[2.0,6.0]
    bandpass = [2.0, 8.0]
    # plotted duration for stream +-det_dur
    # det_dur=25.0
    det_dur = 20.0

    # carica template
    stemp_num = str(template_num)
    # print "template_num == ", stemp_num
    for file in glob.glob('./template/' + stemp_num + '.*'):
        # print 'file :' + file
        st_temp += read(file)
    # print "Stream ...", st_temp
    # npanels=len(st_temp)
    npanels = 27
    tt = Trace()
    # calcola il tempo di riferimento minimo dello stream dei template
    # print st_temp
    refT = min([tt.stats.starttime for tt in st_temp])
    # print "refT == ", refT
    st_temp.sort(keys=["starttime"])[0:30]

    # carica tracce 24h continue
    # for file1 in glob.glob('./24h/' + sday + '.?????' + '.' + channel[0]):
    for file1 in glob.glob('./24h/' + sday + '.????' + '.' + channel[
        0]):
        st_cont += read(file1)
    try:
        for file2 in glob.glob('./24h/' + sday + '.???' + '.' + channel[
            0]):
            if os.path.isfile(file2):
                st_cont += read(file2)
    except OSError:
        pass
    try:
        for file3 in glob.glob('./24h/' + sday + '.?????' + '.' +
                                       channel[0]):
            if os.path.isfile(file3):
                st_cont += read(file3)
    except OSError:
        pass
    d = {}

    # define travel time file for each template (travel time files for synchronizing CFTs are obtained
    # running calcTT01.py
    # travel_file="%s%s.ttimes" % (travel_dir, template_num)

    # store ttimes info in a dictionary
    # with open(travel_file, "r") as ttim:
    # d = dict(x.rstrip().split(None, 1) for x in ttim)
    ttmin = [float(str(x)) for x in list(d.values())]
    # print ttmin
    # min_time=min(ttmin[:])
    # print d
    # print " min_time == ", min_time

    # select filtered template
    # st1_temp = st_temp.select(channel=channel[0])
    st1_temp = st_temp.sort(keys=["starttime"])[0:27].select(channel=channel[0])
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
        sta_lat = sta_coord.get(sstat)[0]
        sta_lon = sta_coord.get(sstat)[1]
        # print "sta_lat, sta_lon === ", eve_lat, eve_lon, eve_dep, sta_lat, sta_lon
        ori = calc_timeshift(eve_lat, eve_lon, eve_dep, sta_lat, sta_lon)
        # print "ori == ", ori

        chan = tt.stats.channel
        netwk = tt.stats.network
        idchan = netwk + "." + sstat + ".." + chan
        idchan_dic = netwk + "." + sstat + "." + chan
        # print "idchan_dic == ", idchan_dic
        # print "d[idchan_dic] == ", d[idchan_dic]
        st1_cont = Stream()
        st1_cont = st_cont.select(id=idchan)
        try:
            if st1_cont.__nonzero__():
                # print(st1_cont.__str__(extended=True))
                st1_cont.merge(method=1, fill_value=0)
                tc = st1_cont[0]
                tc.detrend('constant')
                tc.filter("bandpass", freqmin=bandpass[0], freqmax=bandpass[1], zerophase=True)
                # print "Trace Id.i= ", st1_cont[0]
                # tc.trim(starttime=UTCDateTime(detection_otime), endtime=UTCDateTime(detection_otime)+2*det_dur)
                tc.trim(starttime=UTCDateTime(detection_otime), endtime=UTCDateTime(detection_otime) + 1 * det_dur)
                # print "tc.stats.starttime == ", tc.stats.starttime
                # print "detection_otime == ", UTCDateTime(detection_otime)
                try:
                    tssize = tc.data.size
                    if tssize > 0:
                        ind_tmin = int((ori) / tc.stats.delta) + 1
                        ind_tmax = int((ori + 5.0) / tc.stats.delta) + 1
                        magg = max(abs(tc.data))
                        tm = tc.copy()
                        aa = abs(tc.data[ind_tmin:ind_tmax])
                        # amaxad=max(abs(tc.data[ind_tmin:ind_tmax]))
                        if len(aa) == 0:
                            amaxad = 1
                        else:
                            amaxad = max(aa)
                        # amaxad=max(abs(tc.data))
                        amaxat = max(abs(tt.data))
                        # print "amaxat, amaxad == ", amaxat, amaxad
                        amaxmul = amaxad / amaxat
                        # print "amaxmul == ", amaxmul
                        magd = mag_detect(magt, amaxat, amaxad)
                        smagd = str("%4.2f" % magd)
                        smagd = "Md=" + smagd + " Mt=" + str(magt)
                        mag = amaxad
                        tlen = UTCDateTime(detection_otime).timestamp

                        # per plottare statistiche
                        sta_net = str(tc.stats.network) + "." + str(tc.stats.station)
                        sta_chan = str(tc.stats.channel)
                        indices = [i for i, s in enumerate(lines) if sta_net in s and sta_chan in s]
                        # print indices
                        if indices:
                            cc_sta = lines[int(indices[0])][3]
                            shift_sta = lines[int(indices[0])][4]

                        # tt.trim(starttime=UTCDateTime(template_otime), endtime=UTCDateTime(ttbegin)+5, fill_value=0)
                        ############################################################################################
                        ## Plot detected events versus template
                        ############################################################################################
                        tad = np.arange(0, (tt.stats.npts / tt.stats.sampling_rate), tt.stats.delta)
                        t = np.arange(0, tc.stats.npts / tc.stats.sampling_rate, tc.stats.delta)
                        axarray[count - 1].plot(tad + ori, amaxmul * tt.data, 'r', lw=1.0)
                        axarray[count - 1].plot(t, tc.data, 'k', lw=1.5)
                        # axarray[count-1].plot(tad+ori, amaxmul*tt.data, 'r', lw=1.5)
                        axarray[count - 1].text(det_dur * 0.05, 0.65 * magg, smagd, fontsize=11)
                        axarray[count - 1].text(det_dur * 0.25, 0.65 * magg, 'avcc= ' + str(avcc), fontsize=11)
                        axarray[count - 1].text(det_dur * 0.7, 0.65 * magg, 'thre= ' + str(thre), fontsize=11)
                        # axarray[count-1].text(det_dur, 0.65*magg, 'template_num=' + str(template_num), fontsize=11)
                        axarray[count - 1].text(det_dur * 0.9, 0.65 * magg, tc.stats.station + '.' + tc.stats.channel,
                                                fontsize=11)
                        axarray[count - 1].text(det_dur * 0.4, 0.65 * magg, 'cc= ' + cc_sta, fontsize=11)
                        axarray[count - 1].text(det_dur * 0.6, 0.65 * magg, 'shift= ' + shift_sta, fontsize=11)
                        # axarray[count-1].text(det_dur, -0.9*magg, 'template_time=' + str(tt.stats.starttime), fontsize=11)
                        # axarray[count-1].text(det_dur, -0.9*magg, 'template_time=' + str(tt.stats.starttime-ori), fontsize=11)
                        axarray[count - 1].text(det_dur * 0.05, -0.9 * magg, 'detection_time=' + str(detection_otime),
                                                fontsize=11)
                        if count == 1:
                            axarray[count - 1].text(det_dur * 0.3, 1.4 * magg,
                                                    'Mag= ' + str(dcmag[detection_num]) + '    ''Detection = ' + str(
                                                        detection_num) + '    ' + 'Template = ' + str(template_num),
                                                    fontsize=16)
                        plt.xlabel('Time [s]')
                    else:
                        print("warning no stream is found")
                except OSError:
                    pass
            else:
                print("warning no stream is found")
        except OSError:
            pass
        # plt.ylim(-amaxat,amaxat)
        # fig = plt.gcf()
        # fig.set_size_inches(11.93, 15.98,11.93)
        ################################################################################
        # Set Flag_Save_Figure=1 to save png files or Flag_Save_Figure=0 to show results
        ################################################################################
    if Flag_Save_Figure == 0:
        plt.show()
    if Flag_Save_Figure == 1:
        fig = plt.gcf()
        fig.set_size_inches(11.93, 15.98)
        # fig.set_size_inches(15.98, 11.93)
        # print "sfig===", sfig
        # outfile=sday + sstat + "." + str(detection_num) + "." + str(template_num).zfill(3) + ".png"
        outfile = str(detection_otime) + "." + str(template_num).zfill(3) + ".10Hz.png"
        print("outfile", outfile)

        fig.savefig(outfile, dpi=300, papertype='a3')
        fig.clf()
