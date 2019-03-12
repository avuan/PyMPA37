.. image:: pympa_logo1.png
    :width: 100px
    :align: left
    :alt: pympa_logo1.png
    :target: https://github.com/avuan/PyMPA37/releases

Tutorial 
======================================

This tutorial is designed to give you an overview of the capabilities and
implementation of the PyMPA Python package.

Downloading Seismological Data
------------------------------
To download seismological data from EIDA (European Integrated Data Archive) servers: 
Waveform broad band data from european broad band seismic stations are available 
from many others European Institutions. To download seismological data from EIDA (European Integrated Data Archive) servers
and inventory data in STATIONXML format many examples can be found in ObsPy examples.

PyMPA requires, continuous data and stations inventories.
EIDA servers can easily release continuous data from permanent networks and the corresponding 
inventories. The examples in the subdirectory input.download_data.dir show the python scripts
that allows the download. 

In the case your data come from other sources, PyMPA through ObsPy libraries
is able to manage most of the seismological data formats (MSEED, SAC, SEISAN, SEGY, etc..). 
An inventory data file including station information needs to be created by modifying an existing 
StationXML file and read that with ObsPy to perform the rest.

PyMPA does not use databases and prefers to store single channel daily continuous data in archieves.

(:doc:`download_data </sub/input.download_data>`) download data

Create Templates
----------------
(:doc:`create_template </sub/input.create_templates>`) create templates

A Python script create_templates.py is used to trim templates from continuous data
stored in an archive. Generally, we use S-wave travel times to cut events before and after arrivals.
Take care that a high sampling rate could result in memory consumption
and prolonged times of execution. Input data should be decimated a priori accordingly with your needs.
Check example for running create_templates.py at https://github.com/avuan/PyMPA37/tree/master/input.create_templates.dir 

Needed files:

- Events in a catalog: e.g. templates.zmap (quakeml or zmap format) see ObsPy for format

.. include:: ../../input.create_templates.dir/templates.zmap
   :literal:

- Suitable velocity model for computing travel times

.. include:: ../../input.create_templates.dir/aquila_kato.tvel
   :literal:

- Station inventory (format consistent with ObsPy read_inventory routine see https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.read_inventory.html)

- Days to process: one column file including days to process e.g. lista1

.. include:: ../../input.create_templates.dir/lista1
   :literal:

- Set parameters: e.g. trim.par

.. include:: ../../input.create_templates.dir/trim.par 
   :literal: 
   
- Directory i.e. ./24h where 24h continuous data are stored

- Output dir i.e. ./template (find trimmed time series)

Note that input and output file names, inventories, template catalogs, velocity models are recalled also in the next steps and remain almost fixed. Parameters in files .par could change.
  

Check Template Quality
----------------------
(:doc:`template_check </sub/input.template_check>`) select good templates

Evaluating template quality allows to input only a good signal to noise ratio avoiding artifacts resulting in unwanted detections. The selection is based on Kurtosis method (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kurtosis.html) supposing that the waveform is simmetrically trimmed at the first S-wave arrival. Kurtosis evaluates if the time distrbution of amplitudes is simmetric or not excluding data having a low signal to noise ratio.

Check example running test_kurtosis1.py at https://github.com/avuan/PyMPA37/tree/master/input.template_check.dir
After running kurtosis based selection, templates are separated in two subdirectories "bad"(red waveforms see figure below) and "good" (black waveforms see figure below)

.. image:: ./figure/bad.png
    :width: 400px
    :align: center

.. image:: ./figure/bad4.png
    :width: 400px
    :align: center

.. image:: ./figure/bad7.png
    :width: 400px
    :align: center

.. image:: ./figure/good.png
    :width: 400px
    :align: center

.. image:: ./figure/good1.png
    :width: 400px
    :align: center

.. image:: ./figure/good3.png
    :width: 400px
    :align: center


Calculate Travel Times
----------------------
(:doc:`calculate_ttimes </sub/input.calculate_ttimes>`) calculate travel times

Travel time calculation is based on Java TauP Toolkit as implemented in ObsPy (https://docs.obspy.org/packages/obspy.taup.html) 

Needed files:

- Events in a catalog: e.g. templates.zmap (quakeml or zmap format) see ObsPy for format

- Suitable velocity model for computing travel times

- Station inventory (format consistent with ObsPy read_inventory routine see https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.read_inventory.html)

- Days to process: one column file including days to process e.g. lista1

- Set parameters: e.g. times.par

.. include:: ../../input.calculate_ttimes.dir/times.par
   :literal:

- Input directory i.e. ./template where trimmed templates are found

- Output dir i.e. ./ttimes (find moveout times from different channels used to synchronize cross-correlation functions)

Note that input and output file names, inventories, template catalogs, velocity models are recalled also in the next steps and remain almost fixed. Parameters in files .par could change.



Running PyMPA
-------------
(:doc:`main.pympa_chunks_channel_limit </sub/main.pympa_chunks_channel_limit>`) run pympa

Template matching code, using cross-correlation based on well located eventsr. The code is embarassingly parallel and different templates/days can be run on different cores. We do not provide the scripts to parallelize jobs preferring to leave to the user to find the best strategy to accomplish the task.
Three different versions of the matching code are provdided:

1) main.pympa.py (working on chunks of 24 hours and using all the available channels). This version costs a lot in terms of memory usage and time.
2) main.pympa_channel_limit.py (working on a limited number of channels close to the template epicenter). It saves computation time reducing also the memory consumption.
3) main.pympa_chunks_channel_limit.py (working on daily chunks and with a reduced number of channels). This is the fastest version saving also a lot of memory and avoiding possible memory leakage occurring sometimes in version 1.
 
The input files are the same by using the different versions. The only input file that changes is the parameters24 input file that in the case of version 2 and 3 has added respectively 1 or 2 input lines.

.. include:: ../../main.pympa.dir/parameters24
   :literal:

.. include:: ../../main.pympa_channel_limit.dir/parameters24
   :literal: 

.. include:: ../../main.pympa_chunks_channel_limit.dir/parameters24
   :literal:

 
Needed input files in common between versions 1, 2, 3:

- Events in a catalog: e.g. templates.zmap (quakeml or zmap format) see ObsPy for format
- Suitable velocity model for computing travel times
- Station inventory (format consistent with ObsPy read_inventory routine see https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.read_inventory.html)
- Days to process: one column file including days to process e.g. lista1
- Set parameters: e.g. parameters24.par
- Input directory ./template where trimmed templates are found
- Input directory ./24h where 24 hours continuous waveforms are stores  
- Input directory /ttimes (find moveout times from different channels used to synchronize cross-correlation functions)


Output:

- Output files .cat, .stats, .stats.mag, .except (details on the output )

Detections (.cat)

Detections are listed in a .cat file (e.g. 200.100723.cat)

.. include:: ../../main.pympa_chunks_channel_limit.dir/200.100723.cat
   :literal:

Columns in 200.100723.cat file are:

- Template number corresponding to the python line index in file templates.zmap (event catalog)
- UTC Date and Time (2010-07-23T02:20:29.833321Z) Time precision selection is possible in parameters24 input file
- Magnitude estimated as in Peng and Zhao (2009). The magnitude of the detected event is calculated as the median value of the maximum amplitude ratios for all channels between the template and detected event, assuming that a 10-fold increase in amplitude corresponds to a one-unit increase in magnitude.
- Average cross-correlation estimated from the channels that concurred to the detection. This value is estimated using a time shift between the channels that optimized the stacked CFT.
- Threshold value (ratio between the amplitude of the CFT stacking and the daily MAD Median Absolute Deviation). The higher the threshold the most probable the detection. This value is estimated using a time shift between the channels that optimized the stacked CFT.
- Average cross-correlation estimated from the channels that concurred to the detection at no time shift.
- Threshold value (ratio between the amplitude of the CFT stacking and the daily MAD Median Absolute Deviation). No time shift of the signal cross-correlation functions is allowed.
- Number of channels for which the cross-correlation is over a certain lower bound (e.g. 0.35)

Single Channel Statistics is listed in a .stats file (e.g. 200.100723.stats)

.. include:: ../../main.pympa_chunks_channel_limit.dir/200.100723.stats
   :literal:

Columns in 200.100723.stats file are:

- Network.Station
- Channel
- Cross-correlation value at no time shift 
- Cross-correlation value with time shift (nsamples) as in column 5
- Time shift in nsamples (e.g. -1.0 means that the shift is equal to 0.05 at 20Hz sampling rate)

At the end of each trace id you find other parameters related to the detection in part repeating the detection parameters 
in .cat file and in part related to the cross-correlations values over some limits (0.3 - 0.5 - 0.7 - 0.9).

- date, template_num, detection_num, date&time, template_magnitude, detection_magnitude, threshold_fixed, MAD, ave_crosscc, threshold_record, ave_crosscc_0, threshold_record_0, num_channels_gt0.3, num_channels_gt0.5, num_channels_gt0.7, num_channels_gt0.9
- 100723 201 0 2010-07-23T22:20:57.712239Z 1.51 0.07 9.0 0.0342230009997 0.486 14.193 0.301 8.796 11.0 5.0 2.0 0.0

.. include:: ../../main.pympa_chunks_channel_limit.dir/200.100723.stats.mag
   :literal:

Columns in 200.100723.stats.mag file are:

- Station.Channel Mag.
- date, template_num, detection_num, date&time, template_magnitude, detection_magnitude, threshold_fixed, MAD, ave_crosscc, threshold_record, ave_crosscc_0, threshold_record_0, num_channels_gt0.3, num_channels_gt0.5, num_channels_gt0.7, num_channels_gt0.9
- 100723 201 0 2010-07-23T22:20:57.712239Z 1.51 0.07 9.0 0.0342230009997 0.486 14.193 0.301 8.796 11.0 5.0 2.0 0.0

Note that input and output file names, inventories, template catalogs, velocity models are recalled also in the next steps and remain almost fixed. Parameters in files .par could change.




Output Processing
-----------------
(:doc:`output.process_detections </sub/output.process_detections>`) controls multiple detections in short time windows


Verify Detections
-----------------
(:doc:`output.verify_detection </sub/output.verify_detection>`) visual verification of events 




References
----------

Shelly, D. R., G. C. Beroza, and S. Ide (2007). Non-volcanic tremor and low
frequency earthquake swarms, Nature 446, 305–307.

Peng, Z., and P. Zhao (2009). Migration of early aftershocks following the
2004 Parkfield earthquake, Nature Geosci. 2, 877–881.

Yang, H., L. Zhu, and R. Chu (2009). Fault-plane determination of the
18 April 2008 Mount Carmel, Illinois, earthquake by detecting and
relocating aftershocks, Bull. Seismol. Soc. Am. 99, 3413–3420.

Kato, A., K. Obara, T. Igarashi, H. Tsuruoka, S. Nakagawa, and N. Hirata
(2012). Propagation of slow slip leading up to the 2011 Mw 9.0
Tohoku-Oki earthquake, Science 335, 705–708.

Zhang, M., and L. Wen (2015). An effective method for small event detection:
Match and locate (M&L), Geophys. J. Int. 200, 1523–1537.

Krischer, L., T. Megies, R. Barsch, M. Beyreuther, T. Lecocq, C. Caudron,
and J. Wassermann (2015). ObsPy: A bridge for seismology into the
scientific Python ecosystem, Comput. Sci. Discov. 8, no. 1, 014003,
doi: 10.1088/1749-4699/8/1/014003.
