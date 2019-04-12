Running PyMPA
-------------

The input files are those prepared in advance by using pre-processing tools. The only input file that changes is the parameters24 input file.

.. include:: ../../../main.pympa.dir/parameters24
   :literal:


Needed input files:

- Events in a catalog: e.g. templates.zmap (quakeml or zmap format) see ObsPy for format
- Suitable velocity model for computing travel times
- Station inventory (format consistent with ObsPy read_inventory routine see https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.read_inventory.html)
- Days to process: one column file including days to process e.g. lista1
- Set parameters: e.g. parameters24
- Input directory ./template where trimmed templates are found
- Input directory ./24h where 24 hours continuous waveforms are stores
- Input directory /ttimes (find moveout times from different channels used to synchronize cross-correlation functions)


Output:

- Output files .cat, .stats, .stats.mag, .except (details on the output )

Detections (.cat)

Detections are listed in a .cat file (e.g. 200.100723.cat)

.. include:: ../../../main.pympa.dir/200.100723.cat
   :literal:

Columns in 200.100723.cat file are:

- 1.Template number corresponding to the python line index in file templates.zmap (event catalog)
- 2.UTC Date and Time (2010-07-23T02:20:29.833321Z) Time precision selection is possible in parameters24 input file
- 3.Magnitude estimated as in Peng and Zhao (2009). The magnitude of the detected event is calculated as the median value of the maximum amplitude ratios for all channels between the template and detected event, assuming that a 10-fold increase in amplitude corresponds to a one-unit increase in magnitude.
- 4.Average cross-correlation estimated from the channels that concurred to the detection. This value is estimated using a time shift between the channels that optimized the stacked CFT.
- 5.Threshold value (ratio between the amplitude of the CFT stacking and the daily MAD Median Absolute Deviation). The higher the threshold the most probable the detection. This value is estimated using a time shift between the channels that optimized the stacked CFT.
- 6.Average cross-correlation estimated from the channels that concurred to the detection at no time shift.
- 7.Threshold value (ratio between the amplitude of the CFT stacking and the daily MAD Median Absolute Deviation). No time shift of the signal cross-correlation functions is allowed.
- 8.Number of channels for which the cross-correlation is over a certain lower bound (e.g. 0.35)

Single Channel Statistics is listed in a .stats file (e.g. 200.100723.stats)

.. include:: ../../../main.pympa.dir/200.100723.stats
   :start-line: 0
   :end-line: 26
   :literal:

Columns in 200.100723.stats file are:

- 1.Network.Station
- 2.Channel
- 3.Cross-correlation value at no time shift
- 4.Cross-correlation value with time shift (nsamples) as in column 5
- 5.Time shift in nsamples (e.g. -1.0 means that the shift is equal to 0.05 at 20Hz sampling rate)

At the end of each trace id you find other parameters related to the detection in part repeating the detection parameters
in .cat file and in part related to the cross-correlations values over some limits (0.3 - 0.5 - 0.7 - 0.9).

- date, template_num, detection_num, date&time, template_magnitude, detection_magnitude, threshold_fixed, MAD, ave_crosscc, threshold_record, ave_crosscc_0, threshold_record_0, num_channels_gt0.3, num_channels_gt0.5, num_channels_gt0.7, num_channels_gt0.9
- 100723 201 0 2010-07-23T22:20:57.712239Z 1.51 0.07 9.0 0.0342230009997 0.486 14.193 0.301 8.796 11.0 5.0 2.0 0.0

.. include:: ../../../main.pympa.dir/200.100723.stats.mag
   :start-line: 0
   :end-line: 26
   :literal:

Columns in 200.100723.stats.mag file are:

- 1.Station.Channel Mag.
- 2.date, template_num, detection_num, date&time, template_magnitude, detection_magnitude, threshold_fixed, MAD, ave_crosscc, threshold_record, ave_crosscc_0, threshold_record_0, num_channels_gt0.3, num_channels_gt0.5, num_channels_gt0.7, num_channels_gt0.9
- 3.100723 201 0 2010-07-23T22:20:57.712239Z 1.51 0.07 9.0 0.0342230009997 0.486 14.193 0.301 8.796 11.0 5.0 2.0 0.0

Note that input and output file names, inventories, template catalogs, velocity models are recalled also in the next steps and remain almost fixed. Parameters in files .par could change.

