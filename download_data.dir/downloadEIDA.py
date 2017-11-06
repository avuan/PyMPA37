#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import time

from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.event import read_events


def call_bulk(net, sta, chann, start, stop, chuncklength):
    t1 = (UTCDateTime(start)).timestamp
    t3 = (UTCDateTime(stop)).timestamp
    t2 = t1 + chuncklength

    while (t2<=t3):
        print("station == ", sta)

        for net in networks:

            for chann in channels:
                try:
                    bulk = [(net, sta, "*", chann, UTCDateTime(t1), UTCDateTime(t2))]
                    print("bulk == ", bulk)
                    yy = str(UTCDateTime(t1).year)
                    mm = str(UTCDateTime(t1).month).zfill(2)
                    dd = str(UTCDateTime(t1).day).zfill(2)

                    newfile = inp_dir + yy + mm + dd + "." + sta + "." + chann[0:2]
                    client.get_waveforms_bulk(bulk, filename=newfile)
                    time.sleep(2)
                except Exception:
                    pass
# 
client = Client("INGV")
networks = ["MN", "IV"]
stations = ["AQU", "ARRO", "CAMP", "FEMA", "FDMO", "GIGS", "LNSS", "MMO1",
            "NRCA", "RM33", "TERO", "T1243", "AQT1", "SMA1", "OFFI", "GUMA",
            "CESI", "MOMA", "CESX"]
channels = ["EH*", "HH*", "HN*"]
# 24h as seconds 
chuncklength = 86400

# output directory
inp_dir = "./24h/"

# if ichoice == 0 select a period from start to stop
# if ichoice == 1 select time windows from a templates.zmap catalog
ichoice = 0

if ichoice == 0:
    start = "2012-06-26T00:00:00.000"
    stop = "2012-06-28T00:00:00.000"

    for sta in stations:

        for net in networks:

            for chann in channels:
                call_bulk(net, sta, chann, start, stop, chuncklength)
            start = start + chuncklength
elif ichoice == 1:
    cat = read_events("nu.zmap", format="ZMAP")
    ncat = len(cat)

    for iev in range(0, ncat):
        ot = cat[iev].origins[0].time.datetime
        ot1 = UTCDateTime(ot)
        yy = ot1.year
        mm = ot1.month
        dd = ot1.day
        start = str(yy) + "-" + str(mm) + "-" + str(dd) + "T00:00:00.000"
        stop = UTCDateTime(start) + chuncklength

        for sta in stations:

            for net in networks:

                for chann in channels:
                    call_bulk(net, sta, chann, start, stop, chuncklength)
