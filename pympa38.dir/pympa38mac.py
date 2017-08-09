#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 2016/08/23 Version 34 - parameters24 input file needed
#
# August 2015 - August 2016 Alessandro Vuan
# Code for the detection of microseismicity based on cross correlation of template events.
# The code exploits multiple cores to speed up time (nwork parameter has to be set to use multiple cores)
# 
# Method's references:
# Kato A, Obara K, Igarashi T, Tsuruoka H, Nakagawa S, Hirata N (2012) Propagation
# of slow slip leading up to the 2011 Mw 9.0 Tohoku-Oki earthquake. Science
# doi:10.1126/science.1215141 
#
# Software Requirements: the following dependencies are needed (check import and from statements below)
# Python "obspy" package installed via Anaconda with all numpy and scipy packages
# Python "math" libraries
# Python "bottleneck" utilities to speed up numpy array operations
#
# import useful libraries

import glob
import math as M
import os
import os.path
import sys
from math import log10

import bottleneck as bn
import numpy as np
import obspy
import obspy.signal
from obspy import read
from obspy.core import *
from obspy.core.event import *
from obspy.signal import *
from obspy.signal.trigger import coincidence_trigger
from obspy.signal.util import *


# LIST OF USEFUL FUNCTIONS

def TrimFillOneDay(tc, iiyear, iimonth, iiday, iihour, iimin, iisec):
    starttime_sec = UTCDateTime(iiyear, iimonth, iiday, iihour, iimin, iisec).timestamp
    starttime_sec = starttime_sec - 86400
    endtime_sec = starttime_sec + 86400
    t00 = UTCDateTime(starttime_sec)
    t11 = UTCDateTime(endtime_sec)
    tc.trim(starttime=t00, endtime=t11, pad=True, fill_value=0)
    return tc


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def xcorr(x, y):
    N = len(x)
    M = len(y)
    meany = bn.nanmean(y)
    stdy = bn.nanstd(np.asarray(y))
    tmp = rolling_window(x, M)
    with np.errstate(divide='ignore'):
        c = bn.nansum((y - meany) * (tmp - np.reshape(bn.nanmean(tmp, -1), (N - M + 1, 1))), -1) / (
            M * bn.nanstd(tmp, -1) * stdy)
        c[M * bn.nanstd(tmp, -1) * stdy == 0] = 0
        return c


def processInput(itemp, nn, ss, ich, stream_df):
    global stream_cft
    #   itemp = template number, nn =  network code, ss = station code, ich = channel code, stream_df = Stream() object as defined in obspy library
    temp_file = "%s.%s.%s..%s.mseed" % (str(itemp), nn, ss, ich)
    finpt = "%s%s" % (temp_dir, temp_file)
    if os.path.isfile(finpt):
        tsize = 0
        try:
            tsize = os.path.getsize(finpt)
            if tsize > 0:
                # print "ok template exist and not empty"
                st_temp = Stream()
                st_temp = read(finpt)
                tt = st_temp[0]
                # continuous data are stored in stream_df
                sc = stream_df.select(station=ss, channel=ich)
                if sc.__nonzero__():
                    tc = sc[0]
                    fct = xcorr(tc.data, tt.data)
                    stats = {'network': tc.stats.network, 'station': tc.stats.station, 'location': '',
                             'channel': tc.stats.channel, 'starttime': tc.stats.starttime,
                             'npts': len(fct), 'sampling_rate': tc.stats.sampling_rate,
                             'mseed': {'dataquality': 'D'}}
                    trnew = Trace(data=fct, header=stats)
                    tc = trnew.copy()
                    # TrimFillOneDay(tc,iyear,imonth,iday)
                    stream_cft += Stream(traces=[tc])
                    # tc.write(newfile, format = "MSEED")
                else:
                    print("warning no stream is found")
            else:
                print("warning template event is empty")
        except OSError:
            pass


def quality_cft(trac):
    std_trac = bn.nanstd(abs(trac.data))
    return std_trac


def stack(stall, df, tstart, npts, stdup, stddown):
    std_trac = np.empty(len(stall))
    """
    Function to stack traces in a stream with different trace.id and different starttime
    but the same number of datapoints. Returns a trace having as starttime
    the earliest startime within the stream
    """
    for itr, tr in enumerate(stall):
        std_trac[itr] = quality_cft(tr)
    avestd = bn.nanmean(std_trac[0:])
    avestdup = avestd * stdup
    avestddw = avestd * stddown

    for itr, tr in enumerate(stall):
        if (std_trac[itr] >= avestdup or std_trac[itr] <= avestddw):
            stall.remove(tr)
            print("removed Trace n Stream = ...", tr, std_trac[itr], avestd)
    itr = len(stall)
    tdat = bn.nansum([tr.data for tr in stall], axis=0) / itr
    sta = "STACK"
    cha = "BH"
    net = "XX"
    header = {'network': net, 'station': sta,
              'channel': cha, 'starttime': tstart, 'sampling_rate': df, 'npts': npts}
    tt = Trace(data=tdat, header=header)
    return tt


def check_singlechannelCFT(stall, stCC, trg, tstda, sample_tol, cc_threshold, nch_min, day, itemp, itrig, f1):
    """
    The function check_singlechannelCFT compute the maximum CFT's values
    at each trigger time and counts the number of channels having higher cross-correlation
    nch, cft_ave, cftratio are re-evaluated on the basis of +/- 2 sample approximation.
    Statistics are written in stat files
    """
    # important parameters: a sample_tolerance less than 2 results in wrong magnitude
    sample_tolerance = sample_tol
    single_channelCFT = cc_threshold
    #
    trigger_time = trg['time']
    tcft = stCC[0]
    t0_tcft = tcft.stats.starttime
    trigger_shift = trigger_time.timestamp - t0_tcft.timestamp
    trigger_sample = round(trigger_shift / tcft.stats.delta)
    max_sct = np.empty(len(stall))
    max_trg = np.empty(len(stall))
    std_sct = np.empty(len(stall))
    max_ind = np.empty(len(stall))
    chan_sct = np.chararray(len(stall), 12)
    nch = 0
    for icft, tsc in enumerate(stall):
        # get cft amplitude value at corresponding trigger and store it in
        # check for possible 2 sample shift and eventually change trg['cft_peaks']
        chan_sct[icft] = tsc.stats.network + "." + tsc.stats.station + " " + tsc.stats.channel
        max_sct[icft] = max(tsc.data[trigger_sample - sample_tolerance:trigger_sample + sample_tolerance + 1])
        max_ind[icft] = np.argmax(tsc.data[trigger_sample - sample_tolerance:trigger_sample + sample_tolerance + 1])
        max_ind[icft] = sample_tolerance - max_ind[icft]
        max_trg[icft] = tsc.data[trigger_sample:trigger_sample + 1]
    nch = (max_sct > single_channelCFT).sum()
    if nch >= nch_min:
        nch09 = (max_sct > 0.9).sum()
        nch07 = (max_sct > 0.7).sum()
        nch05 = (max_sct > 0.5).sum()
        nch03 = (max_sct > 0.3).sum()
        #    print("nch ==", nch03, nch05, nch07, nch09)
        cft_ave = bn.nanmean(max_sct[:])
        cftratio = cft_ave / tstda
        cft_ave_trg = bn.nanmean(max_trg[:])
        cftratio_trg = cft_ave_trg / tstda
        max_sct = max_sct.T
        max_trg = max_trg.T
        chan_sct = chan_sct.T
        # str11 = "%s %s %s %s %s %s %s %s %s %s %s %s %s \n" % (day[0:6], str(itemp), str(itrig),
        #    trigger_time, tstda, cft_ave, cftratio, cft_ave_trg, cftratio_trg, nch03, nch05, nch07, nch09)
        # str11 = "%s %s %s %s %s %s %s %s \n" % ( nch03, nch04, nch05, nch06, nch07, nch08, cft_ave, cftratio )
        # f1.write(str11)
        for idchan in range(0, len(max_sct)):
            str22 = "%s %s %s %s \n" % (chan_sct[idchan].decode(), max_trg[idchan], max_sct[idchan], max_ind[idchan])
            f1.write(str22)
    else:
        nch = 1
        cft_ave = 1
        cftratio = 1
        cft_ave_trg = 1
        cftratio_trg = 1
        nch03 = 1
        nch05 = 1
        nch07 = 1
        nch09 = 1

    return nch, cft_ave, cftratio, cft_ave_trg, cftratio_trg, nch03, nch05, nch07, nch09


def mag_detect(magt, amaxt, amaxd):
    """
    mag_detect(mag_temp,amax_temp,amax_detect)
    Returns the magnitude of the new detection by using the template/detection amplitude trace ratio
    and the magnitude of the template event
    """
    amaxr = amaxt / amaxd
    magd = magt - log10(amaxr)
    return magd

def reject_Moutliers(data, m=1.):
    ndata = len(data)
    nonzeroind = np.nonzero(data)[0]
    nzlen = len(nonzeroind)
    print("nonzeroind ==", nonzeroind)
    data = data[nonzeroind]
    print("data ==", data)
    datamed = np.median(data)
    print("datamed ==", datamed)
    d = np.abs(data - datamed)
    mdev = 2 * np.median(d)
    print("d, mdev ==", d, mdev)
    if mdev == 0:
        inds = np.arange(nzlen)
        print("inds ==", inds)
        data[inds] = datamed
    else:
        s = d / mdev
        inds = np.where(s <= m)
        print("inds ==", inds)
    return data[inds]


# read 'parameters24' file to setup useful variables

with open('parameters24') as fp:
    data = fp.read().splitlines()
stations = data[19].split(" ")
channels = data[20].split(" ")
networks = data[21].split(" ")
lowpassf = float(data[22])
highpassf = float(data[23])
sample_tol = int(data[24])
cc_threshold = float(data[25])
nch_min = int(data[26])
temp_length = float(data[27])
UTC_prec = int(data[28])
cont_dir = "./" + data[29] + "/"
temp_dir = "./" + data[30] + "/"
travel_dir = "./" + data[31] + "/"
day_list = str(data[32])
ev_catalog = str(data[33])
start_itemp = int(data[34])
stop_itemp = int(data[35])
factor_thre = int(data[36])
stdup = float(data[37])
stddown = float(data[38])
# ---------------------------------
# set time precision for UTCDATETIME
UTCDateTime.DEFAULT_PRECISION = UTC_prec

# read Catalog of Templates Events
cat = Catalog()
cat = read_events(ev_catalog, format="ZMAP")
ncat = len(cat)

# read template from standard input
# startTemplate = input("INPUT: Enter Starting template ")
# stopTemplate = input("INPUT: Enter Ending template ")
# print("OUTPUT: Running from template", startTemplate,  " to ", stopTemplate)
t_start = start_itemp
t_stop = stop_itemp

# t_start = int(startTemplate)
# t_stop = int(stopTemplate)

fname = day_list

# array of days is built deleting last line character (/newline) ls -1 command include a newline character at the end
with open(fname) as fl:
    days = [line[:-1] for line in fl]
    print(days)
fl.close()

# loop over days
"""
initialise stt as a stream of templates
and stream_df as a stream of continuous waveforms
"""
stt = Stream()
stream_df = Stream()
stream_cft = Stream()
for day in days:
    # settings to cut exactly 24 hours file from without including previous/next day 
    iday = "%s" % (day[4:6])
    imonth = "%s" % (day[2:4])
    print("imonth ==", imonth)
    iyear = "20%s" % (day[0:2])
    iiyear = int(iyear)
    iimonth = int(imonth)
    iiday = int(iday)
    iihour = 23
    iimin = 59
    iisec = 0
    # Prepare stream of continuous waveforms 
    stream_df.clear()
    tc = Trace()
    bandpass = [lowpassf, highpassf]
    finpc = "%s%s.*.???" % (cont_dir, str(day))
    for ic, filenamec in enumerate(glob.glob(finpc)):
        st = read(filenamec)
        st.merge(method=1, fill_value=0)
        tc = st[0]
        stat = tc.stats.station
        chan = tc.stats.channel
        tc.detrend('constant')
        # ensuring that 24h continuous trace starts at 00 hour 00 minut 00.0 seconds
        TrimFillOneDay(tc, iiyear, iimonth, iiday, iihour, iimin, iisec)
        tc.filter("bandpass", freqmin=bandpass[0], freqmax=bandpass[1], zerophase=True)
        #       store detrended and filtered continuous data in a Stream
        #       this is needed to evaluate a new detection magnitude if triggers are found
        #       before applying the convolution verify that template event exists
        stream_df += Stream(traces=[tc])
    for itemp in range(t_start, t_stop):
        stt.clear()
        fout = "%s.%s.cat" % (str(itemp), day[0:6])
        f = open(fout, 'w+')
        print("itemp == ...", str(itemp))
        ot = cat[itemp].origins[0].time
        mt = cat[itemp].magnitudes[0].mag
        lon = cat[itemp].origins[0].longitude
        lat = cat[itemp].origins[0].latitude
        dep = cat[itemp].origins[0].depth
        # amplitude info on templates important to detect magnitude of new events
        inplist = "%s%s.??.*..???.mseed" % (temp_dir, str(itemp))

        for filename in glob.glob(inplist):
            csize = os.path.getsize(filename)
            if csize > 0:
                # store in Stream template events
                stt += read(filename)
        # check number of zeroes within the template and avoid using it in case perc < 75%
        for tr in stt:
            npts = tr.stats.npts
            ck0 = np.count_nonzero(tr.data)
            perc = 100 * float(ck0) / float(npts)
            if perc <= 75:
                stt.remove(tr)
        ntl = len(stt)
        amaxat = np.empty(ntl)
        # for each template event
        # md=np.empty(ntl)
        md = np.zeros(ntl)
        damaxat = {}
        # reference time to be used for retrieving time synchronization
        refT = min([tr.stats.starttime for tr in stt])
        for il, tr in enumerate(stt):
            amaxat[il] = max(abs(tr.data))
            sta_t = tr.stats.station
            cha_t = tr.stats.channel
            tid_t = "%s.%s" % (sta_t, cha_t)
            damaxat[tid_t] = amaxat[il]
        d = {}
        # define travel time file for each template (travel time files for synchronizing CFTs are obtained
        # running calcTT01.py
        travel_file = "%s%s.ttimes" % (travel_dir, str(itemp))
        print("travel_file = ", travel_file)
        # store ttimes info in a dictionary
        with open(travel_file, "r") as ttim:
            d = dict(x.rstrip().split(None, 1) for x in ttim)
            ttim.close()
        # print(d)
        # find minimum time to recover origin time
        time_values = [float(v) for v in d.values()]
        min_time_value = min(time_values)
        print("min_time_value == ", min_time_value)
        min_time_key = [k for k, v in d.items() if v == str(min_time_value)]
        print("key, mintime == ", min_time_key, min_time_value)

        stream_cft.clear()
        for nn in networks:
            for ss in stations:
                for ich in channels:
                    processInput(itemp, nn, ss, ich, stream_df)

        stnew = Stream()
        stall = Stream()
        tr = Trace()
        tt = Trace()
        tc_cft = Trace()
        tsnew = UTCDateTime()

        # seconds in 24 hours
        h24 = 86400
        nfile = len(stream_cft)
        Tstart = np.empty(nfile)
        Tend = np.empty(nfile)
        tdif = np.empty(nfile)
        for idx, tc_cft in enumerate(stream_cft):
            # get station name from trace
            # s=tc_cft.id
            sta = tc_cft.stats.station
            chan = tc_cft.stats.channel
            net = tc_cft.stats.network
            delta = tc_cft.stats.delta

            npts = h24 / delta
            s = "%s.%s.%s" % (net, sta, chan)
            # read S-wave travel time for synchro  (network.station.channel = column[0] arr.s correction = column[1] )
            # Scol=columns[1]
            # print(d)
            tdif[idx] = float(d[s])
        tdifmin = min(tdif[0:])
        for idx, tc_cft in enumerate(stream_cft):
            # get stream starttime
            Tstart[idx] = tc_cft.stats.starttime + tdif[idx]
            # waveforms should have the same number of npts
            # and should be synchronized to the S-wave travel time
            secs = h24 + 60
            Tend[idx] = Tstart[idx] + secs
            check_npts = (Tend[idx] - Tstart[idx]) / tc_cft.stats.delta
            ts = UTCDateTime(Tstart[idx], precision=UTC_prec)
            te = UTCDateTime(Tend[idx], precision=UTC_prec)
            stall += tc_cft.trim(starttime=ts, endtime=te, nearest_sample=True, pad=True, fill_value=0)

        tstart = min([tr.stats.starttime for tr in stall])
        df = stall[0].stats.sampling_rate
        npts = stall[0].stats.npts

        # compute mean cross correlation from the stack of CFTs (see stack function)
        ccmad = stack(stall, df, tstart, npts, stdup, stddown)

        # compute standard deviation of abs(ccmad)
        # tstda=np.std(abs(ccmad.data))

        ccm = ccmad.data[ccmad.data != 0]
        tstda = bn.nanmedian(abs(ccm))

        # define threshold as 9 times std  and quality index
        thresholdD = (factor_thre * tstda)
        # Trace ccmad is stored in a Stream
        stCC = Stream(traces=[ccmad])
        # Run coincidence trigger on a single CC trace resulting from the CFTs stack

        # essential threshold parameters
        # Cross correlation thresholds
        xcor_cut = thresholdD
        thr_on = thresholdD
        thr_off = thresholdD - 0.15 * thresholdD
        thr_coincidence_sum = 1.0
        similarity_thresholds = {"BH": thr_on}
        trigger_type = None
        triglist = coincidence_trigger(trigger_type, thr_on, thr_off, stCC, thr_coincidence_sum, trace_ids=None,
                                       similarity_thresholds=similarity_thresholds, delete_long_trigger=False,
                                       trigger_off_extension=3.0, details=True)
        ntrig = len(triglist)

        tt = np.empty(ntrig)
        cs = np.empty(ntrig)
        nch = np.empty(ntrig)
        cft_ave = np.empty(ntrig)
        cftratio = np.empty(ntrig)
        cft_ave_trg = np.empty(ntrig)
        cftratio_trg = np.empty(ntrig)
        nch3 = np.empty(ntrig)
        nch5 = np.empty(ntrig)
        nch7 = np.empty(ntrig)
        nch9 = np.empty(ntrig)
        mm = np.empty(ntrig)
        timex = UTCDateTime()
        fout1 = "%s.%s.stats" % (str(itemp), day[0:6])
        f1 = open(fout1, 'w+')
        for itrig, trg in enumerate(triglist):
            if tdifmin == min_time_value:
                tt[itrig] = trg['time'] + tdifmin
            elif tdifmin != min_time_value:
                diff_time = min_time_value - tdifmin
                tt[itrig] = trg['time'] + diff_time
            cs[itrig] = trg['coincidence_sum']
            cft_ave[itrig] = trg['cft_peak_wmean']
            cftratio[itrig] = trg['cft_peaks'][0] / tstda
            traceID = trg['trace_ids']
            # check single channel CFT
            nch[itrig], cft_ave[itrig], cftratio[itrig], cft_ave_trg[itrig], cftratio_trg[itrig], nch3[itrig], nch5[
                itrig], nch7[itrig], nch9[itrig] = check_singlechannelCFT(stall, stCC, trg, tstda, sample_tol,
                                                                          cc_threshold, nch_min, day, itemp, itrig, f1)
            if int(nch[itrig]) >= nch_min:
                nn = len(stream_df)
                # nn=len(stt)
                amaxac = np.zeros(nn)
                md = np.zeros(nn)

                # for each trigger, detrended, and filtered continuous data channels are trimmed and amplitude useful to estimate 
                # magnitude is measured.        
                damaxac = {}
                timex = UTCDateTime(tt[itrig])
                for il, tc in enumerate(stream_df):
                    ss = tc.stats.station
                    ich = tc.stats.channel
                    netwk = tc.stats.network
                    if stt.select(station=ss, channel=ich).__nonzero__():
                        ttt = stt.select(station=ss, channel=ich)[0]
                        s = "%s.%s.%s" % (netwk, ss, ich)
                        # print " s ==", s
                        if tdifmin < 0:
                            timestart = timex + abs(tdifmin) + (
                                UTCDateTime(ttt.stats.starttime).timestamp - UTCDateTime(refT).timestamp)
                        elif tdifmin > 0:
                            timestart = timex - abs(tdifmin) + (
                                UTCDateTime(ttt.stats.starttime).timestamp - UTCDateTime(refT).timestamp)

                        timend = timestart + temp_length
                        ta = Trace()
                        ta = tc.copy()
                        ta.trim(starttime=timestart, endtime=timend, pad=True, fill_value=0)
                        amaxac[il] = max(abs(ta.data))
                        tid_c = "%s.%s" % (ss, ich)
                        damaxac[tid_c] = amaxac[il]
                        if float(damaxac[tid_c]) != 0 and float(damaxat[tid_c]) != 0:
                            print " float(damaxat[tid_c], float(damaxac[tid_c] == ", float(damaxat[tid_c]), float(damaxac[tid_c])
                            md[il] = mag_detect(mt, float(damaxat[tid_c]), float(damaxac[tid_c]))
                mdr = reject_Moutliers(md, 1)
                mm[itrig] = round(np.mean(mdr), 2)
                cft_ave[itrig] = round(cft_ave[itrig], 3)
                cftratio[itrig] = round(cftratio[itrig], 3)
                cft_ave_trg[itrig] = round(cft_ave_trg[itrig], 3)
                cftratio_trg[itrig] = round(cftratio_trg[itrig], 3)
                str33 = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % (
                    day[0:6], str(itemp), str(itrig), str(UTCDateTime(tt[itrig])), str(mm[itrig]), str(mt),
                    str(nch[itrig]),
                    str(tstda), str(cft_ave[itrig]), str(cftratio[itrig]), str(cft_ave_trg[itrig]),
                    str(cftratio_trg[itrig]), str(nch3[itrig]), str(nch5[itrig]), str(nch7[itrig]), str(nch9[itrig]))
                f1.write(str33)
                str1 = "%s %s %s %s %s %s %s %s\n" % (
                    str(itemp), str(UTCDateTime(tt[itrig])), str(mm[itrig]), str(cft_ave[itrig]), str(cftratio[itrig]),
                    str(cft_ave_trg[itrig]), str(cftratio_trg[itrig]), str(int(nch[itrig])))
                f.write(str1)
    f1.close()
f.close()
