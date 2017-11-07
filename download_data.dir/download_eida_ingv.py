#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import time
import os
from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime


client = Client("INGV")
networks = ["MN", "IV"]
#
# stations=["AQU", "ARRO", "CAMP", "FEMA", "FDMO", "GIGS",
# "LNSS", "MMO1", "NRCA", "RM33", "TERO", "T1243", "AQT1", "SMA1",
# "OFFI", "GUMA", "CESI", "MOMA", "CESX"]

stations = ["AQU"]
channels = ["EHZ", "EHN", "EHE", "HHZ", "HHN", "HHE", "HNZ", "HNN", "HNE"]

start = "2016-01-01T00:00:00.000"
stop = "2016-01-02T00:00:00.000"

# 24h as seconds
chuncklength = 86400

# output directory
inp_dir = "./24h/"

# output directory
if not os.path.exists(inp_dir):
    os.makedirs(inp_dir)
# ---------------------------do not change below

for sta in stations:
    t1 = (UTCDateTime(start)).timestamp
    t3 = (UTCDateTime(stop)).timestamp
    t2 = t1 + chuncklength

    while t2 <= t3:
        print("station == ", sta)

        for net in networks:

            for chann in channels:

                try:
                    bulk = [(net, sta, "*", chann, UTCDateTime(
                        t1), UTCDateTime(t2))]
                    print("bulk == ", bulk)
                    yy = str(UTCDateTime(t1).year)
                    mm = str(UTCDateTime(t1).month).zfill(2)
                    dd = str(UTCDateTime(t1).day).zfill(2)

                    newfile = inp_dir + yy + mm + dd +\
                        "." + sta + "." + chann[0:3]
                    client.get_waveforms_bulk(bulk, filename=newfile)
                    time.sleep(2)
                except Exception:
                    pass

        t1 = t1 + chuncklength
        t2 = t2 + chuncklength
