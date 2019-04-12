Download data
*************

Downloading data from EIDA servers is easily performed by using ObsPy tools
and routines. Please see the scripts at https://github.com/avuan/PyMPA37/blob/master/input.download_data.dir/

download_data.py and download_inventory.py allows downloading continuous data and inventories for the selected stations
and time period.

For some networks or stations please check before data availability or if there are some restrictions and a token
or userid is needed. Some EIDA servers when receiving from the same user many requests in a short time
close the door to the user. That is the reason why we introduced a pause command in our scripts, this avoids
the shutdown of the connection.

Here below a simple code to download continuous data:

.. code-block:: python

    #!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import time
import os
from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime

client = Client("INGV")
networks = ["IV"]
#
stations=["MURB", "NARO"]

channels = ["EH*", "HH*"]

start = "2012-02-05T00:00:00.000"
stop = "2012-02-06T00:00:00.000"

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
                    time.sleep(1)
                    pass

        t1 = t1 + chuncklength
        t2 = t2 + chuncklength

