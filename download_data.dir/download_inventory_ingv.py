#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from obspy.clients.fdsn import Client
from obspy import Catalog, UTCDateTime
from obspy import read
import time

client = Client("INGV")
#networks=["MN", "IV"]

starttime = UTCDateTime("2009-01-01")
endtime = UTCDateTime("2017-01-02")

inventory = client.get_stations(
    starttime=starttime, endtime=endtime,
    network="AI", sta="*", channel="*", filename="inv.ingv.ai", format='xml') 
print(inventory)
