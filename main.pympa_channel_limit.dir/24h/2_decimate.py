#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
import glob
import json
import obspy.signal
import numpy as np
import matplotlib.pyplot as plt

from math import log10
from obspy.core import *
from obspy.signal import *
from obspy.signal.util import *
from obspy.imaging import *
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import *
from obspy.taup.taup_create import *
from obspy.io.mseed.util import get_timing_and_data_quality
from obspy import Catalog, UTCDateTime
from obspy.core.event import * 
from obspy.io.zmap import *


# setup station list
#bandpass=[2.0,5.0]
stations=["ATBU", "ATCA", "ATFO", "ATLO", "ATMI", "ATPC", "ATPI", "ATRE", "ATSC", "ATVA", "ATVO", "BADI", "CDCA", "MURB", "PIEI"]

stations=["APEC", "NARO", , , , "FRON", , ,
          , "SSFR", ,  "FOSV", "ATCC", "PE3", "PARC",
          


"BADI", "ATBU", "MURB", "ATLO", "ATPC", "ATVO", "ATFO", "ATPI", "ATSC"
networks=["IV"]
channels=["HH", "EH"]
components=["E", "N", "Z"]
inp_dir="/Volumes/Taboo/test1/original/"
#out_dir="/Volumes/Taboo/test1/24h/"
nlen=len(inp_dir)
for sta in stations:
    inp_files=inp_dir + "2013*." + sta + ".??"
    for file in glob.glob(inp_files):
        print(file)
        st=Stream()
        st=read(file)
        sta=st[0].stats.station
        print("sta ==", sta)
        for comp in components:
            st1=Stream()
            st1=st.select(component=comp)
            #print(st1)
            #chan=st1[0].stats.channel
            if len(st1) > 0:
                chan=st1[0].stats.channel
                st1.merge(method = 1, fill_value = 0)
                st1.detrend('constant')
                samp = st1[0].stats.sampling_rate
                if (samp/25.).is_integer():
                    dec_factor=int(samp/25)
                    st1.decimate(dec_factor, strict_length=False, no_filter=True)
                    newsamp = st1[0].stats.sampling_rate
                    newfile= file[0:21]+ "24h/" + file[nlen+2:nlen+8]+ "." + sta + "." + chan 
                    print(newfile)
                    if st1.__nonzero__():
                       st1[0].data=st1[0].data.astype(np.int32)
                       st1.write(newfile, format = "MSEED", encoding='STEIM1')
        #try:
            #os.remove(file)
        #except OSError:
            #pass

