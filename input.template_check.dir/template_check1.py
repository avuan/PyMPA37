from obspy import *
from scipy.stats import kurtosis
from subprocess import call
import matplotlib.pyplot as plt
import os.path
import numpy as np
from obspy.signal.filter import envelope

# FlagPrint==0 no printout, FlagPrint==1 check waveforms
FlagPrint = 1

st = Stream()
tr = Trace()
call(["rm", "-r", "./bad"])
call(["rm", "-r", "./good"])
call(["mkdir", "./bad"])
call(["mkdir", "./good"])
channels = ["EHE", "EHN", "EHZ", "HHE", "HHN", "HHZ"]
stations = ["NRCA"]
networks = ["IV", "MN"]

for s in stations:

    for c in channels:

        for n in networks:

            for i in range(0, 570):
                inpfile = str(i) + "." + n + "." + s + ".." + c + ".mseed"
                try:
                    if os.path.isfile(inpfile):
                        st = read(inpfile)
                        tr = st[0]
                        dt = tr.stats.delta
                        npts = tr.stats.npts
                        # print(npts)
                        t = np.linspace(0, (npts - 1) * dt, npts)
                        # print(t)
                        amaxa = max(tr.data)

                        # first test
                        # check if we are in the seismogram coda
                        # maximum absolute amplitude is at the beginning
                        index_max = np.argmax(abs(tr.data))
                        pp = int(npts * 0.3)

                        if index_max <= pp:
                            indm = 0
                        elif index_max > pp:
                            indm = 1

                        # second test
                        # check if maximum amplitude is in the last 5% of coda points
                        pp95 = int(npts * 0.95)
                        if index_max >= pp95:
                            indm95 = 0
                        elif index_max < pp95:
                            indm95 = 1

                        # tre = envelope(tr.data)
                        k = kurtosis(
                            tr.data,
                            axis=0,
                            fisher=True,
                            bias=True,
                            nan_policy="propagate",
                        )
                        kstring = (
                            "kurt(Pearson) = "
                            + str(k)[0:6]
                            + " ind1 = "
                            + str(indm)
                            + " ind2 = "
                            + str(indm95)
                        )
                        # print(inpfile, k, indm, indm95)

                        if k <= 0.8 or indm == 0 or indm95 == 0:
                            call(["cp", inpfile, "./bad"])
                            print(inpfile, k, indm, indm95)
                        elif k > 0.8 and indm == 1 and indm95 == 1:
                            call(["cp", inpfile, "./good"])

                        if FlagPrint == 1:
                            fig = plt.figure(figsize=(8, 3))
                            ax = fig.add_subplot(111)
                            ax.set_title(inpfile)

                            if k <= 0.8 or indm == 0 or indm95 == 0:
                                ax.plot(t, tr.data, "r")
                            elif k > 0.8 and indm == 1 and indm95 == 1:
                                ax.plot(t, tr.data, "k")
                            ax.set_xlabel("Time[s]")
                            ax.text(0.5, 0.9 * amaxa, kstring, fontsize=15)
                            plt.show()
                except:
                    print("warning file not found ", inpfile)
                    pass
