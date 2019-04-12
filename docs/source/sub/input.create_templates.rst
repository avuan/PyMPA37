Create templates
****************

Actually templates are created by trimming a fixed time window focused on S-wave theoretical travel times. For details on the travel time calculations see file:///Users/vuan/PycharmProjects/PyMPA37/docs/build/html/tutorial.html#calculate-travel-times.

The length of the time window is established by the inter-statiion network distance and the frequency range used.
The user should carefully check to exclude signal deriving from numerical artifacts (e.g. filtering applied to zero padding
time windows having no data), or pre and coda signals not connected with the seismic perturbation investigated (e.g. LFEs, earthquakes, icequakes etc...)
In the next versions the trimming will allow for selecting variable length P and S-waves.

Needed files:

- Events in a catalog: e.g. templates.zmap (quakeml or zmap format) see ObsPy for more details. An example of zmap file format is given here below (ZMAP is a simple 10 column CSV file (technically TSV) format for basic catalog data. It originates from ZMAP, a Matlab® based earthquake statistics package (see [Wiemer2001]). 
- Columns represent: longitude, latitude, year, month, day, magnitude, depth, hour, minute, second - note that fields have to be separated by tabs.

.. include:: ../../../input.create_templates.dir/zmap.txt
   :start-line: 0
   :end-line: 10
   :literal:

- Suitable velocity model for computing travel times

.. include:: ../../../input.create_templates.dir/aquila_kato.tvel
   :start-line: 0
   :end-line: 20
   :literal:

- Station inventory (format consistent with ObsPy read_inventory routine see https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.read_inventory.html)

- Days to process: one column file including days to process e.g. lista1

.. include:: ../../../input.create_templates.dir/lista1
   :literal:

- Set parameters: e.g. trim.par

.. include:: ../../../input.create_templates.dir/trim.par
   :literal:
  
- Directory i.e. ./24h where 24h continuous data are stored

- Output dir i.e. ./template (find trimmed time series)

Note that input and output file names, inventories, template catalogs, velocity models are recalled also in the next steps and remain almost fixed. Parameters in files .par could change.

References
Wiemer S. (2001), A software package to analyzeseismicity: ZMAP, Seismol Res Lett 92, 373-382
