#!/usr/bin/env python
#
# import useful libraries
import numpy as np
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.event import read_events
from obspy.io.quakeml.core import _is_quakeml
from obspy.io.zmap.core import _is_zmap


def empty(value):
    try:
        value = float(value)
    except ValueError:
        pass
    return bool(value)


# read 'filter.par' file to setup useful variables

with open("filter.par") as fp:
    data = fp.read().splitlines()

FlagDateTime = int(data[8])
Flag_cc = int(data[9])
inp_file = str(data[10])
utc_prec = int(data[11])
window_length = float(data[12])
min_threshold = float(data[13])
min_nch = int(data[14])
ev_catalog = str(data[15])

print("FlagDateTime == ", FlagDateTime)
print("Flag_cc == ", Flag_cc)
print("Input File == ", inp_file)
print("utc_prec == ", utc_prec)
print("Window_Length == ", window_length)
print("min_threshold == ", min_threshold)
print("min_nch == ", min_nch)
print("ev_catalog == ", ev_catalog)

fp.close()

# read event coordinates from catalog

print("event catalog should be ZMAP or QUAKEML")

if _is_zmap(ev_catalog):
    print("reading ZMAP catalog")

elif _is_quakeml(ev_catalog):
    print("reading QUAKEML catalog")

else:
    print("warning error in reading ZMAP or QUAKEML")

cat = read_events(ev_catalog)

# read 'filter.par' file to setup useful variables
# ppp

# if FlagDateTime==0 output time is in unix seconds 123845678,23
# to be used for plotTM.gmt script
# if FlagDateTime==1 output time is a timestamp e.g 2009-03-30T16:34:25.12
# FlagDateTime=1
# if Flag_cc==0 read cc_ave and cc_sum at trigger time elif Flag_cc ==1
# read the best cc_ave and cc_sum
# within the sample tolerance interval
# Flag_cc=1
# defines directorY of cat.out
# cont_dir=["./"]
# set time precision for UTCDATETIME
# UTCDateTime.DEFAULT_PRECISION=6
# set itime interval before and after detection time to apply for a search
# tt_int=6.0
# min_threshold
# min_nch

# loop over detections
# inp1=input()
fileinp = inp_file
# count input file's number of lines
num_lines = sum(1 for line in open(fileinp))
print(" detections (all) == ", num_lines)
# define empty array for each field data
template_num = np.empty(num_lines)
times = np.empty(num_lines)
mt = np.empty(num_lines)
cc_ave = np.empty(num_lines)
cc_sum = np.empty(num_lines)
cc_ave0 = np.empty(num_lines)
cc_sum0 = np.empty(num_lines)
nch = np.empty(num_lines)

# read daily catalog as an output from running all templates on a
# one day coninuous data
#
f = open(fileinp, "r")

for ind, line in enumerate(f):
    line = line.strip()
    # print(repr(line))
    columns = line.split()
    #    print("ind, columns == ", ind, len(columns))
    template_num[ind] = columns[0]
    times[ind] = float(UTCDateTime(columns[1]))
    #    print("TIMES ==", times[ind])
    mt[ind] = columns[2]

    if Flag_cc == 1:
        cc_ave[ind] = columns[3]
        cc_sum[ind] = columns[4]
        cc_ave0[ind] = columns[5]
        cc_sum0[ind] = columns[6]
    elif Flag_cc == 0:
        cc_ave[ind] = columns[5]
        cc_sum[ind] = columns[6]

    nch[ind] = columns[7]

# search is performed on the number of detections
# define an empty numpy array
#
det = np.chararray(num_lines, 80)
# for each new detection evaluate the time span and define indeces
# corresponding to events that are detected within this time interval

for it, tt in enumerate(times):
    ttmin = tt - window_length
    ttmax = tt + window_length
    # print("ttmin, ttmax ===", ttmin, ttmax)
    indeces = np.where(np.logical_and(times >= ttmin, times <= ttmax))

    # printout indeces
    # print("it, indeces == ", it, indeces)

    # find the maximum ccmad and the corresponding
    # "index" between indeces of events that are within a time interval
    # for each new detection we search for the index of the detection
    # within the same time window that shows the maximum CCMAD

    threshold = np.max(cc_sum[indeces])
    ave_th = np.max(cc_ave[indeces])

    # if multiple detections with the same value for cc_sum are found
    ind_winners = np.argwhere(cc_sum[indeces] == np.max(cc_sum[indeces]))
    nind = len(ind_winners)

    # if multiple detections with cc_ave == 1.0
    if nind > 1 and ave_th == 1.0:
        winner = 0

        # find index for which the cc_sum at sample tolerance 0 is
        # maximum
        for iind in ind_winners.ravel():
            iind1 = np.asarray(indeces).flatten()[iind]

            if cc_sum0[iind1] > winner:
                winner = cc_sum0[iind1]
                index = iind

    else:
        index = np.argmax(cc_sum[indeces])

    # store the event value (0-15390) in the variable det[it]
    # print("it, threshold, index,  == ", it, threshold, indeces[0][index],
    # template_num[indeces[0][index]])

    # print("Flag 1111", template_num[indeces[0][index]],
    # round(times[indeces[0][index]], 2),  mt[indeces[0][index]],
    # cc_ave[indeces[0][index]], cc_sum[indeces[0][index]])
    strtemp = str(template_num[indeces[0][index]])
    ttdet = str(UTCDateTime(times[indeces[0][index]]))
    # print("ttdet == ", ttdet)
    strmt = str(mt[indeces[0][index]])
    strcc_ave = str(cc_ave[indeces[0][index]])
    strcc_sum = str(cc_sum[indeces[0][index]])
    strnch = str(nch[indeces[0][index]])
    # print("it", it)
    seq = (strtemp, ttdet, strmt, strcc_ave, strcc_sum, strnch)
    dets = " ".join(seq)
    # print("dets == ", dets)
    det[it] = dets
    # print("det == ", det[it].decode())
    # find unique detections - avoid duplications - and printout
#
detections = np.unique(det.decode())
nu = len(detections)
print(" detections are reduced to == ", nu)
outfile = open(fileinp + "f1", "w+")

for iu in detections:
    # print(iu)
    outfile.write(iu + "\n")

outfile.close()

# reorder events by time
#
# loop over detections
fileinp1 = fileinp + "f1"
print("fileinp1 == ", fileinp1)
# count input file's number of lines
num_lines1 = sum(1 for line in open(fileinp1))
print("num_lines1 == ", num_lines1)
# define empty array for each field data
t_num = np.empty(num_lines1)
t_tim = np.empty(num_lines1)
t_mag = np.empty(num_lines1)
t_ave = np.empty(num_lines1)
t_sum = np.empty(num_lines1)
t_nch = np.empty(num_lines1)

# read daily catalog as an output from running all templates
# on a one day continuous data
#
f1 = open(fileinp1, "r")

for ind1, line in enumerate(f1):
    line = line.strip()
    # print(repr(line))
    columns1 = line.split()
    # print(columns1[:])
    t_num[ind1] = columns1[0]
    t_tim[ind1] = float(UTCDateTime(columns1[1]))
    t_mag[ind1] = columns1[2]
    t_ave[ind1] = columns1[3]
    t_sum[ind1] = columns1[4]
    t_nch[ind1] = columns1[5]

f1.close()
# sort arrays
t_numr = np.empty(num_lines1)
t_timr = np.empty(num_lines1)
t_magr = np.empty(num_lines1)
t_aver = np.empty(num_lines1)
t_sumr = np.empty(num_lines1)
t_nchr = np.empty(num_lines1)
aa = np.empty(num_lines1)
ttmini = np.empty(num_lines1)
ttmaxi = np.empty(num_lines1)

ireorder = np.argsort(t_tim)
# print ireorder

for irr, jj in enumerate(ireorder):
    # print( "irr, jj == ", irr, jj)
    t_numr[irr] = t_num[jj]
    t_timr[irr] = t_tim[jj]
    # print("t_timr ===", t_timr[irr])
    t_magr[irr] = t_mag[jj]
    t_aver[irr] = t_ave[jj]
    t_sumr[irr] = t_sum[jj]
    t_nchr[irr] = t_nch[jj]

outfile1 = open(fileinp1 + "f2", "w+")

for iu, ti in enumerate(t_timr):
    ttmini[iu] = ti - window_length
    ttmaxi[iu] = ti + window_length
    # print("ti - ttmaxi[iu-1] == ", ti-ttmaxi[iu-1])

    if iu == 0:
        # print("ttmini, ttmaxi ===", UTCDateTime(ttmini[iu]),
        # UTCDateTime(ttmaxi[iu]))
        # find events in a 5 seconds window for each ttmin-ttmax window
        indeces1 = np.where(np.logical_and(t_timr >= ttmini[iu], t_timr <= ttmaxi[iu]))

        # find the maximum ccmad and the corresponding
        # "index" between indeces of events that are within a time interval
        # for each new detection we search for the index of the detection
        # within the same time window that shows the maximum CCMAD
        threshold = np.max(t_sumr[indeces1[0]])
        index1 = np.argmax(t_sumr[indeces1[0]])
        ie = indeces1[0][index1]
        # print("iu, indeces1, index1 == ", iu, indeces1[0], ie)

        ty = UTCDateTime(t_timr[ie]).year
        tmm = UTCDateTime(t_timr[ie]).month
        td = UTCDateTime(t_timr[ie]).day
        th = UTCDateTime(t_timr[ie]).hour
        tminute = UTCDateTime(t_timr[ie]).minute
        tsecond = UTCDateTime(t_timr[ie]).second
        # precision = 10 ** UTCDateTime.DEFAULT_PRECISION
        precision = 10 ** utc_prec
        tmicro = UTCDateTime(t_timr[ie]).microsecond / precision
        tmicro = UTCDateTime(t_timr[ie]).microsecond
        print(
            "tsecond, tmicro == ", tsecond, tmicro, UTCDateTime(t_timr[ie]).microsecond
        )
        correctTime = UTCDateTime(int(ty), int(tmm), int(td)).timestamp
        stringtime = (
            str(ty)
            + " "
            + str(tmm)
            + " "
            + str(td)
            + " "
            + str(th)
            + " "
            + str(tminute)
            + " "
            + str(tsecond)
            + "."
            + str(tmicro).zfill(UTCDateTime.DEFAULT_PRECISION)
        )
        # stringtime=str(ty) + " " + str(tmm) + " " +
        # str(td) + " " + str(th) + " " + str(tminute) + " "
        # + str(tsecond) + "." + str(tmicro)

        if t_sumr[ie] > min_threshold and t_nchr[ie] >= min_nch:

            if FlagDateTime == 1:
                lon = cat[int(t_numr[ie])].origins[0].longitude
                lat = cat[int(t_numr[ie])].origins[0].latitude
                dep = cat[int(t_numr[ie])].origins[0].depth / 1000
                sf = (
                    stringtime
                    + " "
                    + str(t_magr[ie])
                    + " "
                    + str(t_aver[ie])
                    + " "
                    + str(t_sumr[ie])
                    + " "
                    + str(int(t_numr[ie]))
                    + " "
                    + str(lat)
                    + " "
                    + str(lon)
                    + " "
                    + str(dep)
                    + " "
                    + str(t_nchr[ie])
                )
            else:
                t_timr[ie] = t_timr[ie] - correctTime
                sf = (
                    str(t_timr[ie])
                    + " "
                    + str(t_magr[ie])
                    + " "
                    + str(t_sumr[ie])
                    + " "
                    + str(t_nchr[ie])
                )

            outfile1.write(sf + "\n")

    elif iu > 0 and ti > ttmaxi[iu - 1]:
        # print("ttmini, ttmaxi ===", UTCDateTime(ttmini[iu]),
        # UTCDateTime(ttmaxi[iu]))
        # find events in a 5 seconds window for each ttmin-ttmax window
        indeces1 = np.where(np.logical_and(t_timr >= ttmini[iu], t_timr <= ttmaxi[iu]))

        # find the maximum ccmad and the corresponding
        # "index" between indeces of events that are within a time interval
        # for each new detection we search for the index of the detection
        # within the same time window that shows the maximum CCMAD
        threshold = np.max(t_sumr[indeces1[0]])
        index1 = np.argmax(t_sumr[indeces1[0]])
        ie = indeces1[0][index1]
        # print("iu, indeces1, index1 == ", iu, indeces1[0], ie)

        ty = UTCDateTime(t_timr[ie]).year
        tmm = UTCDateTime(t_timr[ie]).month
        td = UTCDateTime(t_timr[ie]).day
        th = UTCDateTime(t_timr[ie]).hour
        tminute = UTCDateTime(t_timr[ie]).minute
        tsecond = UTCDateTime(t_timr[ie]).second
        tmicro = UTCDateTime(t_timr[ie]).microsecond / precision
        tmicro = UTCDateTime(t_timr[ie]).microsecond
        print(
            "tsecond, tmicro == ", tsecond, tmicro, UTCDateTime(t_timr[ie]).microsecond
        )
        correctTime = UTCDateTime(int(ty), int(tmm), int(td)).timestamp
        stringtime = (
            str(ty)
            + " "
            + str(tmm)
            + " "
            + str(td)
            + " "
            + str(th)
            + " "
            + str(tminute)
            + " "
            + str(tsecond)
            + "."
            + str(tmicro).zfill(UTCDateTime.DEFAULT_PRECISION)
        )
        # stringtime=str(ty) + " " + str(tmm) + " " + str(td) + " " +
        # str(th) + " " + str(tminute) + " " + str(tsecond) + "." +
        # str(tmicro)

        if t_sumr[ie] > min_threshold and t_nchr[ie] >= min_nch:

            if FlagDateTime == 1:
                lon = cat[int(t_numr[ie])].origins[0].longitude
                lat = cat[int(t_numr[ie])].origins[0].latitude
                dep = cat[int(t_numr[ie])].origins[0].depth / 1000
                sf = (
                    stringtime
                    + " "
                    + str(t_magr[ie])
                    + " "
                    + str(t_aver[ie])
                    + " "
                    + str(t_sumr[ie])
                    + " "
                    + str(int(t_numr[ie]))
                    + " "
                    + str(lat)
                    + " "
                    + str(lon)
                    + " "
                    + str(dep)
                    + " "
                    + str(t_nchr[ie])
                )
            else:
                t_timr[ie] = t_timr[ie] - correctTime
                sf = (
                    str(t_timr[ie])
                    + " "
                    + str(t_magr[ie])
                    + " "
                    + str(t_sumr[ie])
                    + " "
                    + str(t_nchr[ie])
                )

            outfile1.write(sf + "\n")

    elif iu > 0 and ti < ttmaxi[iu - 1]:
        print(
            "iur<0 ttmini, ttmaxi >===",
            UTCDateTime(ttmini[iu]),
            UTCDateTime(ttmaxi[iu]),
        )

outfile1.close()
