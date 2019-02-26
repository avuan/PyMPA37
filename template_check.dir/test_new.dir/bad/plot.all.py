from obspy import *
import numpy as np
import matplotlib.pyplot as plt

st = read('*NRCA..HHE.mseed')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for tc in st:
    t = np.arange(0, tc.stats.npts / tc.stats.sampling_rate,
                          tc.stats.delta)
    amaxa = max(abs(tc.data))
    ax.plot(t, tc.data/amaxa, "b-")
plt.show()

