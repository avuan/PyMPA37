#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 2016/08/23 Version 34 - parameters24 input file needed
# 2017/10/27 Version 39 - Reformatted PEP8 Code
# 2017/11/05 Version 40 - Corrections to tdifmin, tstda calculations
# 2019/10/15 Version pympa - xcorr substitued with correlate_template from obspy

# First Version August 2014 - Last October 2017 (author: Alessandro Vuan)

# Code for the detection of microseismicity based on cross correlation
# of template events. The code exploits multiple cores to speed up time
#
# Method's references:
# The code is developed and maintained at
# Istituto Nazionale di Oceanografia e Geofisica di Trieste (OGS)
# and was inspired by collaborating with Aitaro Kato and collegues at ERI.
# Kato A, Obara K, Igarashi T, Tsuruoka H, Nakagawa S, Hirata N (2012)
# Propagation of slow slip leading up to the 2011 Mw 9.0 Tohoku-Oki
# earthquake. Science doi:10.1126/science.1215141
#
# For questions comments and suggestions please send an email to avuan@inogs.it
# The kernel function xcorr used from Austin Holland is modified in pympa
# Recommended the use of Obspy v. 1.2.0 with the substitution of xcorr function with
# correlate_template

# Software Requirements: the following dependencies are needed (check import
# and from statements below)
# Python "obspy" package installed via Anaconda with all numpy and scipy
# packages
# Python "math" libraries
# Python "bottleneck" utilities to speed up numpy array operations
#
# import useful libraries

import os
import os.path
import datetime
from math import log10
from time import perf_counter

import bottleneck as bn
import numpy as np
import pandas as pd
from obspy import read, Stream, Trace
from obspy.core import UTCDateTime
from obspy.core.event import read_events
from obspy.signal.trigger import coincidence_trigger


# from obspy.signal.cross_correlation import correlate_template


# LIST OF USEFUL FUNCTIONS


def listdays(year, month, day, period):
    # create a list of days for scanning by templates
    datelist = pd.date_range(datetime.datetime(year, month, day), periods=period).tolist()
    a = list(map(pd.Timestamp.to_pydatetime, datelist))
    days = []
    for i in a:
        days.append(i.strftime("%y%m%d"))
    return days


def read_parameters(par):
    # read 'parameters24' file to setup useful variables

    with open(par) as fp:
        data = fp.read().splitlines()

        stations = data[23].split(" ")
        print(stations)
        channels = data[24].split(" ")
        print(channels)
        networks = data[25].split(" ")
        print(networks)
        lowpassf = float(data[26])
        highpassf = float(data[27])
        sample_tol = int(data[28])
        cc_threshold = float(data[29])
        nch_min = int(data[30])
        temp_length = float(data[31])
        utc_prec = int(data[32])
        cont_dir = "./" + data[33] + "/"
        temp_dir = "./" + data[34] + "/"
        travel_dir = "./" + data[35] + "/"
        dateperiod = data[36].split(" ")
        ev_catalog = str(data[37])
        start_itemp = int(data[38])
        print("starting template = ", start_itemp)
        stop_itemp = int(data[39])
        print("ending template = ", stop_itemp)
        factor_thre = int(data[40])
        stdup = float(data[41])
        stddown = float(data[42])
        chan_max = int(data[43])
        nchunk = int(data[44])
        return (
            stations,
            channels,
            networks,
            lowpassf,
            highpassf,
            sample_tol,
            cc_threshold,
            nch_min,
            temp_length,
            utc_prec,
            cont_dir,
            temp_dir,
            travel_dir,
            dateperiod,
            ev_catalog,
            start_itemp,
            stop_itemp,
            factor_thre,
            stdup,
            stddown,
            chan_max,
            nchunk,
        )


def trim_fill(tc, t1, t2):
    tc.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    return tc


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def xcorr(x, y):
    n = len(x)
    m = len(y)
    meany = np.nanmean(y)
    stdy = np.nanstd(np.asarray(y))
    tmp = rolling_window(x, m)
    with np.errstate(divide="ignore"):
        c = bn.nansum(
            (y - meany) * (tmp - np.reshape(bn.nanmean(tmp, -1), (n - m + 1, 1))), -1
        ) / (m * bn.nanstd(tmp, -1) * stdy)
        c[m * bn.nanstd(tmp, -1) * stdy == 0] = 0
        return c


def process_input(itemp, nn, ss, ich, stream_df):
    st_cft = Stream()
    # itemp = template number, nn =  network code, ss = station code,
    # ich = channel code, stream_df = Stream() object as defined in obspy
    # library
    temp_file = "%s.%s.%s..%s.mseed" % (str(itemp), nn, ss, ich)
    finpt = "%s%s" % (temp_dir, temp_file)
    if os.path.isfile(finpt):
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
                    # fct = correlate_template(tc.data, tt.data)
                    stats = {
                        "network": tc.stats.network,
                        "station": tc.stats.station,
                        "location": "",
                        "channel": tc.stats.channel,
                        "starttime": tc.stats.starttime,
                        "npts": len(fct),
                        "sampling_rate": tc.stats.sampling_rate,
                        "mseed": {"dataquality": "D"},
                    }
                    trnew = Trace(data=fct, header=stats)
                    tc = trnew.copy()
                    st_cft = Stream(traces=[tc])
                else:
                    print("warning no stream is found")
            else:
                print("warning template event is empty")
        except OSError:
            pass
    return st_cft


def quality_cft(trac):
    std_trac = np.nanstd(abs(trac.data))
    return std_trac


def stack(stall, df, tstart, npts, stdup, stddown, nch_min):
    std_trac = np.empty(len(stall))
    td = np.empty(len(stall))
    """
    Function to stack traces in a stream with different trace.id and
    different starttime but the same number of datapoints.
    Returns a trace having as starttime
    the earliest startime within the stream
    """
    for itr, tr in enumerate(stall):
        std_trac[itr] = quality_cft(tr)
    avestd = np.nanmean(std_trac[0:])
    avestdup = avestd * stdup
    avestddw = avestd * stddown

    for jtr, tr in enumerate(stall):

        if std_trac[jtr] >= avestdup or std_trac[jtr] <= avestddw:
            stall.remove(tr)
            print("removed Trace n Stream = ...", tr, std_trac[jtr], avestd)
            td[jtr] = 99.99
            # print(td[jtr])
        else:
            sta = tr.stats.station
            chan = tr.stats.channel
            net = tr.stats.network
            s = "%s.%s.%s" % (net, sta, chan)
            td[jtr] = float(d[s])
            # print(td[jtr])

    itr = len(stall)
    print("itr == ", itr)
    if itr >= nch_min:
        tdifmin = min(td)
        tdat = np.nansum([tr.data for tr in stall], axis=0) / itr
        sta = "STACK"
        cha = "BH"
        net = "XX"
        header = {
            "network": net,
            "station": sta,
            "channel": cha,
            "starttime": tstart,
            "sampling_rate": df,
            "npts": npts,
        }
        tt = Trace(data=tdat, header=header)

    else:
        tdifmin = None
        sta = "STACK"
        cha = "BH"
        net = "XX"
        header = {
            "network": net,
            "station": sta,
            "channel": cha,
            "starttime": tstart,
            "sampling_rate": df,
            "npts": npts,
        }
        tt = Trace(data=np.zeros(npts), header=header)

    return tt, tdifmin


def csc(
    stall, stcc, trg, tstda, sample_tol, cc_threshold, nch_min, day, itemp, itrig, f1
):
    """
    The function check_singlechannelcft compute the maximum CFT's
    values at each trigger time and counts the number of channels
    having higher cross-correlation
    nch, cft_ave, crt are re-evaluated on the basis of
    +/- 2 sample approximation. Statistics are written in stat files
    """
    # important parameters: a sample_tolerance less than 2 results often
    # in wrong magnitudes
    sample_tolerance = sample_tol
    single_channelcft = cc_threshold
    #
    trigger_time = trg["time"]
    tcft = stcc[0]
    t0_tcft = tcft.stats.starttime
    trigger_shift = trigger_time.timestamp - t0_tcft.timestamp
    trigger_sample = int(round(trigger_shift / tcft.stats.delta))
    max_sct = np.empty(len(stall))
    max_trg = np.empty(len(stall))
    max_ind = np.empty(len(stall))
    chan_sct = np.chararray(len(stall), 12)
    nch = 0

    for icft, tsc in enumerate(stall):
        # get cft amplitude value at corresponding trigger and store it in
        # check for possible 2 sample shift and eventually change
        # trg['cft_peaks']
        chan_sct[icft] = (
            tsc.stats.network + "." + tsc.stats.station + " " + tsc.stats.channel
        )
        tmp0 = trigger_sample - sample_tolerance

        if tmp0 < 0:
            tmp0 = 0
        tmp1 = trigger_sample + sample_tolerance + 1
        max_sct[icft] = max(tsc.data[tmp0:tmp1])
        max_ind[icft] = np.nanargmax(tsc.data[tmp0:tmp1])
        max_ind[icft] = sample_tolerance - max_ind[icft]
        max_trg[icft] = tsc.data[trigger_sample : trigger_sample + 1]
    nch = (max_sct > single_channelcft).sum()

    if nch >= nch_min:
        nch09 = (max_sct > 0.9).sum()
        nch07 = (max_sct > 0.7).sum()
        nch05 = (max_sct > 0.5).sum()
        nch03 = (max_sct > 0.3).sum()
        # print("nch ==", nch03, nch05, nch07, nch09)
        cft_ave = np.nanmean(max_sct[:])
        crt = cft_ave / tstda
        cft_ave_trg = np.nanmean(max_trg[:])
        crt_trg = cft_ave_trg / tstda
        max_sct = max_sct.T
        max_trg = max_trg.T
        chan_sct = chan_sct.T
        # str11 = "%s %s %s %s %s %s %s %s %s %s %s %s %s \n" %
        # (day[0:6], str(itemp), str(itrig),
        # trigger_time, tstda, cft_ave, crt, cft_ave_trg,
        # crt_trg, nch03, nch05, nch07, nch09)
        # str11 = "%s %s %s %s %s %s %s %s \n" % ( nch03, nch04, nch05,
        # nch06, nch07, nch08, cft_ave, crt )
        # f1.write(str11)

        for idchan in range(0, len(max_sct)):
            str22 = "%s %s %s %s \n" % (
                chan_sct[idchan].decode(),
                max_trg[idchan],
                max_sct[idchan],
                max_ind[idchan],
            )
            f1.write(str22)

    else:
        nch = 1
        cft_ave = 1
        crt = 1
        cft_ave_trg = 1
        crt_trg = 1
        nch03 = 1
        nch05 = 1
        nch07 = 1
        nch09 = 1

    return nch, cft_ave, crt, cft_ave_trg, crt_trg, nch03, nch05, nch07, nch09


def mag_detect(magt, amaxt, amaxd):
    """
    mag_detect(mag_temp,amax_temp,amax_detect)
    Returns the magnitude of the new detection by using the template/detection
    amplitude trace ratio
    and the magnitude of the template event
    """
    amaxr = amaxt / amaxd
    magd = magt - log10(amaxr)
    return magd


def reject_moutliers(data, m=1.0):
    nonzeroind = np.nonzero(data)[0]
    nzlen = len(nonzeroind)
    # print("nonzeroind ==", nonzeroind)
    data = data[nonzeroind]
    # print("data ==", data)
    datamed = np.nanmedian(data)
    # print("datamed ==", datamed)
    d = np.abs(data - datamed)
    mdev = 2 * np.median(d)
    # print("d, mdev ==", d, mdev)
    if mdev == 0:
        inds = np.arange(nzlen)
        # print("inds ==", inds)
        data[inds] = datamed
    else:
        s = d / mdev
        inds = np.where(s <= m)
        # print("inds ==", inds)
    return data[inds]


def mad(dmad):
    # calculate daily median absolute deviation
    ccm = dmad[dmad != 0]
    med_val = np.nanmedian(ccm)
    tstda = np.nansum(abs(ccm - med_val) / len(ccm))
    return tstda


start_time = perf_counter()
# read 'parameters24' file to setup useful variables

[
    stations,
    channels,
    networks,
    lowpassf,
    highpassf,
    sample_tol,
    cc_threshold,
    nch_min,
    temp_length,
    utc_prec,
    cont_dir,
    temp_dir,
    travel_dir,
    dateperiod,
    ev_catalog,
    start_itemp,
    stop_itemp,
    factor_thre,
    stdup,
    stddown,
    chan_max,
    nchunk,
] = read_parameters("parameters24")

# set time precision for UTCDATETIME
UTCDateTime.DEFAULT_PRECISION = utc_prec

# read Catalog of Templates Events

cat = read_events(ev_catalog, format="ZMAP")
ncat = len(cat)

# read template from standard input
# startTemplate = input("INPUT: Enter Starting template ")
# stopTemplate = input("INPUT: Enter Ending template ")
# print("OUTPUT: Running from template", startTemplate,  " to ", stopTemplate)
t_start = start_itemp
t_stop = stop_itemp

# loop over days

# generate list of days "
year = int(dateperiod[0])
month = int(dateperiod[1])
day = int(dateperiod[2])
period = int(dateperiod[3])
days = listdays(year, month, day, period)

"""
initialise stt as a stream of templates
and stream_df as a stream of continuous waveforms
"""
stt = Stream()
stream_df = Stream()
stream_cft = Stream()
stall = Stream()
ccmad = Trace()

for day in days:
    # settings to cut exactly 24 hours file from without including
    # previous/next day
    iday = "%s" % (day[4:6])
    imonth = "%s" % (day[2:4])
    print("imonth ==", imonth)
    iyear = "20%s" % (day[0:2])
    iiyear = int(iyear)
    print(iyear, imonth, iday)
    iimonth = int(imonth)
    iiday = int(iday)
    iihour = 23
    iimin = 59
    iisec = 0

    for itemp in range(t_start, t_stop):
        stt.clear()
        # open file containing detections
        fout = "%s.%s.cat" % (str(itemp), day[0:6])
        f = open(fout, "w+")
        print("itemp == ...", str(itemp))
        # open statistics file for each detection
        fout1 = "%s.%s.stats" % (str(itemp), day[0:6])
        f1 = open(fout1, "w+")
        # open file including magnitude information
        fout2 = "%s.%s.stats.mag" % (str(itemp), day[0:6])
        f2 = open(fout2, "w+")
        # open file listing exceptions
        fout3 = "%s.%s.except" % (str(itemp), day[0:6])
        f3 = open(fout3, "w+")

        ot = cat[itemp].origins[0].time
        mt = cat[itemp].magnitudes[0].mag
        lon = cat[itemp].origins[0].longitude
        lat = cat[itemp].origins[0].latitude
        dep = cat[itemp].origins[0].depth

        # read ttimes, select the num_ttimes (parameters,
        # last line) channels
        # and read only these templates
        travel_file = "%s%s.ttimes" % (travel_dir, str(itemp))

        with open(travel_file, "r") as ttim:
            d = dict(x.rstrip().split(None, 1) for x in ttim)
            ttim.close()
            s = d.items()
            v = sorted(s, key=lambda x: (float(x[1])))[0:chan_max]

        vv = [x[0] for x in v]

        for vvc in vv:
            n_net = vvc.split(".")[0]
            n_sta = vvc.split(".")[1]
            n_chn = vvc.split(".")[2]
            filename = "%s%s.%s.%s..%s.mseed" % (
                temp_dir,
                str(itemp),
                str(n_net),
                str(n_sta),
                str(n_chn),
            )
            print(filename)
            stt += read(filename)

        if len(stt) >= nch_min:

            tc = Trace()
            bandpass = [lowpassf, highpassf]

            chunks = []
            h24 = 86400
            chunk_start = UTCDateTime(iiyear, iimonth, iiday)
            end_time = chunk_start + h24

            while chunk_start < end_time:
                chunk_end = chunk_start + h24 / nchunk

                if chunk_end > end_time:
                    chunk_end = end_time
                chunks.append((chunk_start, chunk_end))
                chunk_start += h24 / nchunk

            for t1, t2 in chunks:
                print(t1, t2)
                stream_df.clear()
                for tr in stt:
                    finpc1 = "%s%s.%s.%s" % (
                        cont_dir,
                        str(day),
                        str(tr.stats.station),
                        str(tr.stats.channel),
                    )

                    if os.path.exists(finpc1) and os.path.getsize(finpc1) > 0:

                        try:
                            st = read(finpc1, starttime=t1, endtime=t2)

                            if len(st) != 0:
                                st.merge(method=1, fill_value=0)
                                tc = st[0]
                                stat = tc.stats.station
                                chan = tc.stats.channel
                                tc.detrend("constant")
                                # 24h continuous trace starts 00 h 00 m 00.0s
                                trim_fill(tc, t1, t2)
                                tc.filter(
                                    "bandpass",
                                    freqmin=bandpass[0],
                                    freqmax=bandpass[1],
                                    zerophase=True,
                                )
                                # store detrended and filtered continuous data
                                # in a Stream
                                stream_df += Stream(traces=[tc])

                        except IOError:
                            pass

                if len(stream_df) >= nch_min:

                    ntl = len(stt)
                    amaxat = np.empty(ntl)
                    # for each template event
                    # md=np.empty(ntl)
                    md = np.zeros(ntl)
                    damaxat = {}
                    # reference time to be used for
                    # retrieving time synchronization
                    reft = min([tr.stats.starttime for tr in stt])

                    for il, tr in enumerate(stt):
                        amaxat[il] = max(abs(tr.data))
                        sta_t = tr.stats.station
                        cha_t = tr.stats.channel
                        tid_t = "%s.%s" % (sta_t, cha_t)
                        damaxat[tid_t] = float(amaxat[il])

                    # define travel time file for each template
                    # for synchronizing CFTs are obtained
                    # running calcTT01.py
                    travel_file = "%s%s.ttimes" % (travel_dir, str(itemp))
                    # print("travel_file = ", travel_file)
                    # store ttimes info in a dictionary

                    with open(travel_file, "r") as ttim:
                        d = dict(x.rstrip().split(None, 1) for x in ttim)
                        ttim.close()

                    # print(d)
                    # find minimum time to recover origin time
                    time_values = [float(v) for v in d.values()]
                    min_time_value = min(time_values)
                    # print("min_time_value == ", min_time_value)
                    min_time_key = [k for k, v in d.items() if v == str(min_time_value)]
                    # print("key, mintime == ", min_time_key, min_time_value)

                    # clear global_variable
                    stream_cft.clear()
                    stcc = Stream()

                    for nn in networks:

                        for ss in stations:

                            for ich in channels:
                                stream_cft += process_input(
                                    itemp, nn, ss, ich, stream_df
                                )

                    stall.clear()
                    stcc.clear()
                    stnew = Stream()
                    tr = Trace()

                    tc_cft = Trace()
                    tsnew = UTCDateTime()

                    # seconds in 24 hours

                    nfile = len(stream_cft)

                    tstart = np.empty(nfile)
                    tend = np.empty(nfile)
                    tdif = np.empty(nfile)

                    for idx, tc_cft in enumerate(stream_cft):
                        # get station name from trace
                        sta = tc_cft.stats.station
                        chan = tc_cft.stats.channel
                        net = tc_cft.stats.network
                        delta = tc_cft.stats.delta

                        npts = (h24 / nchunk) / delta
                        s = "%s.%s.%s" % (net, sta, chan)
                        tdif[idx] = float(d[s])

                    for idx, tc_cft in enumerate(stream_cft):
                        # get stream starttime
                        tstart[idx] = tc_cft.stats.starttime + tdif[idx]
                        # waveforms should have the same
                        # number of npts
                        # and should be synchronized to the
                        # S-wave travel time
                        secs = (h24 / nchunk) + 60
                        tend[idx] = tstart[idx] + secs
                        check_npts = (tend[idx] - tstart[idx]) / tc_cft.stats.delta
                        ts = UTCDateTime(tstart[idx], precision=utc_prec)
                        te = UTCDateTime(tend[idx], precision=utc_prec)
                        stall += tc_cft.trim(
                            starttime=ts,
                            endtime=te,
                            nearest_sample=True,
                            pad=True,
                            fill_value=0,
                        )
                    tstart = min([tr.stats.starttime for tr in stall])
                    df = stall[0].stats.sampling_rate
                    npts = stall[0].stats.npts

                    # compute mean cross correlation from the stack of
                    # CFTs (see stack function)

                    ccmad, tdifmin = stack(
                        stall, df, tstart, npts, stdup, stddown, nch_min
                    )
                    print("tdifmin == ", tdifmin)

                    if tdifmin is not None:
                        # compute mean absolute deviation of abs(ccmad)
                        tstda = mad(ccmad.data)

                        # define threshold as 9 times std  and quality index
                        thresholdd = factor_thre * tstda

                        # Trace ccmad is stored in a Stream
                        stcc = Stream(traces=[ccmad])

                        # Run coincidence trigger on a single CC trace
                        # resulting from the CFTs stack
                        # essential threshold parameters
                        # Cross correlation thresholds
                        xcor_cut = thresholdd
                        thr_on = thresholdd
                        thr_off = thresholdd - 0.15 * thresholdd
                        thr_coincidence_sum = 1.0
                        similarity_thresholds = {"BH": thr_on}
                        trigger_type = None
                        triglist = coincidence_trigger(
                            trigger_type,
                            thr_on,
                            thr_off,
                            stcc,
                            thr_coincidence_sum,
                            trace_ids=None,
                            similarity_thresholds=similarity_thresholds,
                            delete_long_trigger=False,
                            trigger_off_extension=3.0,
                            details=True,
                        )
                        ntrig = len(triglist)

                        tt = np.empty(ntrig)
                        cs = np.empty(ntrig)
                        nch = np.empty(ntrig)
                        cft_ave = np.empty(ntrig)
                        crt = np.empty(ntrig)
                        cft_ave_trg = np.empty(ntrig)
                        crt_trg = np.empty(ntrig)
                        nch3 = np.empty(ntrig)
                        nch5 = np.empty(ntrig)
                        nch7 = np.empty(ntrig)
                        nch9 = np.empty(ntrig)
                        mm = np.empty(ntrig)
                        timex = UTCDateTime()
                        tdifmin = min(tdif[0:])

                        for itrig, trg in enumerate(triglist):
                            # tdifmin is computed for contributing channels
                            # within the stack function
                            #

                            if tdifmin == min_time_value:
                                tt[itrig] = trg["time"] + min_time_value

                            elif tdifmin != min_time_value:
                                diff_time = min_time_value - tdifmin
                                tt[itrig] = trg["time"] + diff_time + min_time_value

                            cs[itrig] = trg["coincidence_sum"]
                            cft_ave[itrig] = trg["cft_peak_wmean"]
                            crt[itrig] = trg["cft_peaks"][0] / tstda
                            # traceID = trg['trace_ids']
                            # check single channel CFT
                            [
                                nch[itrig],
                                cft_ave[itrig],
                                crt[itrig],
                                cft_ave_trg[itrig],
                                crt_trg[itrig],
                                nch3[itrig],
                                nch5[itrig],
                                nch7[itrig],
                                nch9[itrig],
                            ] = csc(
                                stall,
                                stcc,
                                trg,
                                tstda,
                                sample_tol,
                                cc_threshold,
                                nch_min,
                                day,
                                itemp,
                                itrig,
                                f1,
                            )

                            if int(nch[itrig]) >= nch_min:
                                nn = len(stream_df)
                                # nn=len(stt)
                                amaxac = np.zeros(nn)
                                md = np.zeros(nn)

                                # for each trigger, detrended,
                                # and filtered continuous
                                # data channels are trimmed and
                                # amplitude useful to
                                # estimate magnitude is measured.
                                damaxac = {}
                                mchan = {}
                                timestart = UTCDateTime()
                                timex = UTCDateTime(tt[itrig])

                                for il, tc in enumerate(stream_df):
                                    ss = tc.stats.station
                                    ich = tc.stats.channel
                                    netwk = tc.stats.network

                                    if stt.select(
                                        station=ss, channel=ich
                                    ).__nonzero__():
                                        ttt = stt.select(station=ss, channel=ich)[0]
                                        s = "%s.%s.%s" % (netwk, ss, ich)
                                        # print " s ==", s
                                        uts = UTCDateTime(ttt.stats.starttime).timestamp
                                        utr = UTCDateTime(reft).timestamp
                                        if tdifmin <= 0:
                                            timestart = (
                                                timex + abs(tdifmin) + (uts - utr)
                                            )

                                        elif tdifmin > 0:
                                            timestart = (
                                                timex - abs(tdifmin) + (uts - utr)
                                            )

                                        timend = timestart + temp_length
                                        ta = tc.copy()
                                        ta.trim(
                                            starttime=timestart,
                                            endtime=timend,
                                            pad=True,
                                            fill_value=0,
                                        )
                                        amaxac[il] = max(abs(ta.data))
                                        tid_c = "%s.%s" % (ss, ich)
                                        damaxac[tid_c] = float(amaxac[il])
                                        dct = damaxac[tid_c]
                                        dtt = damaxat[tid_c]
                                        if dct != 0 and dtt != 0:
                                            md[il] = mag_detect(
                                                mt, damaxat[tid_c], damaxac[tid_c]
                                            )
                                            mchan[tid_c] = md[il]
                                            str00 = "%s %s\n" % (tid_c, mchan[tid_c])
                                            f2.write(str00)

                                mdr = reject_moutliers(md, 1)
                                mm[itrig] = round(np.mean(mdr), 2)
                                cft_ave[itrig] = round(cft_ave[itrig], 3)
                                crt[itrig] = round(crt[itrig], 3)
                                cft_ave_trg[itrig] = round(cft_ave_trg[itrig], 3)
                                crt_trg[itrig] = round(crt_trg[itrig], 3)
                                str33 = (
                                    "%s %s %s %s %s %s %s %s %s "
                                    "%s %s %s %s %s %s %s\n"
                                    % (
                                        day[0:6],
                                        str(itemp),
                                        str(itrig),
                                        str(UTCDateTime(tt[itrig])),
                                        str(mm[itrig]),
                                        str(mt),
                                        str(nch[itrig]),
                                        str(tstda),
                                        str(cft_ave[itrig]),
                                        str(crt[itrig]),
                                        str(cft_ave_trg[itrig]),
                                        str(crt_trg[itrig]),
                                        str(nch3[itrig]),
                                        str(nch5[itrig]),
                                        str(nch7[itrig]),
                                        str(nch9[itrig]),
                                    )
                                )
                                f1.write(str33)
                                f2.write(str33)
                                str1 = "%s %s %s %s %s %s %s %s\n" % (
                                    str(itemp),
                                    str(UTCDateTime(tt[itrig])),
                                    str(mm[itrig]),
                                    str(cft_ave[itrig]),
                                    str(crt[itrig]),
                                    str(cft_ave_trg[itrig]),
                                    str(crt_trg[itrig]),
                                    str(int(nch[itrig])),
                                )
                                f.write(str1)
                    else:
                        str_except2 = "%s %s %s %s %s\n" % (
                            day[0:6],
                            str(itemp),
                            str(t1),
                            str(t2),
                            " num. correlograms lower than nch_min",
                        )
                        f3.write(str_except2)
                        pass
                else:
                    str_except1 = "%s %s %s %s %s\n" % (
                        day[0:6],
                        str(itemp),
                        str(t1),
                        str(t2),
                        " num.  24h channels lower than nch_min",
                    )
                    f3.write(str_except1)
                    pass
        else:
            str_except0 = "%s %s %s\n" % (
                day[0:6],
                str(itemp),
                " num.  templates lower than nch_min",
            )
            f3.write(str_except0)
            pass

        f1.close()
        f2.close()
        f3.close()
        f.close()
print(" elapsed time ", perf_counter() - start_time, " seconds")
