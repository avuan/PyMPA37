from obspy import *
from scipy.stats import kurtosis
from subprocess import call
import matplotlib.pyplot as plt 
import numpy as np
from obspy.signal.filter import envelope
# FlagPrint==0 no printout, FlagPrint==1 check waveforms 
FlagPrint=1
st = Stream()
tr = Trace()
call(['rm', '-r', './bad'])
call(['rm', '-r', './good'])
call(['mkdir', './bad'])
call(['mkdir', './good'])
channels=['HHE','HHN', 'HHZ']
stations=['NRCA']
for s in stations:
    for c in channels: 
        for i in range(0, 574):
            inpfile = str(i) + '.IV.' + s + '..' + c + '.mseed'
            st = read(inpfile)
            tr = st[0]
            dt = tr.stats.delta
            npts = tr.stats.npts
            print(npts)
            t = np.linspace(0, (npts-1)*dt, npts)
            print(t)
            amaxa = max(tr.data)
            # check if we are in the seismogram coda 
            # maximum absolute amplitude is at the beginninga
            index_max = np.argmax(abs(tr.data))
            pp = int(npts * 0.3) 
            if index_max<=pp:
                indm = 0
            elif index_max>pp:
                indm = 1
            #tre = envelope(tr.data)
            k = kurtosis(tr.data, axis=0, fisher=True, bias=True, nan_policy='propagate')
            kstring = 'kurtosis (Pearson) = ' + str(k) + ' coda ind = ' + str(indm)
            print(inpfile, k, indm)
            if k<=0.8 or indm==0:
                call(['cp', inpfile, './bad'])
            elif k>0.8 and indm==1:
                call(['cp', inpfile, './good'])
            if FlagPrint==1:
                fig = plt.figure(figsize=(8,3))
                ax = fig.add_subplot(111)
                ax.set_title(inpfile)
                if k<=0.8 or indm==0:
                    ax.plot(t, tr.data, 'r')
                elif k>0.8 and indm==1:
                    ax.plot(t, tr.data, 'k')
                ax.set_xlabel('Time[s]')
                ax.text(.5, 0.9*amaxa, kstring, fontsize=15)
                plt.show()
