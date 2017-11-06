#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime


client = Client("INGV")

starttime = UTCDateTime("2001-01-01")
endtime = UTCDateTime("2017-06-30")

inventory = client.get_stations(
    starttime=starttime, endtime=endtime,
    network="IV,MN", sta="*", channel="*", filename="inv.ingv", format='xml') 
print(inventory)
