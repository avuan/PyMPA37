Create templates
****************

Actually templates are created by trimming a fixed time window focused on S-wave theoretical travel times. 
The length of the time window is established by the inter-statiion network distance and the frequency range used.
The user should carefully check to exclude signal deriving from numerical artifacts (e.g. filtering applied to zero padding
time windows having no data), or pre and coda signals not connected with the seismic perturbation investigated (e.g. LFEs, earthquakes, icequakes etc...)
In the next versions the trimming will allow for selecting variable length P and S-waves.

Needed files:

- Events in a catalog: e.g. templates.zmap (quakeml or zmap format) see ObsPy for format

.. include:: ../../../input.create_templates.dir/templates.zmap
   :literal:

- Suitable velocity model for computing travel times

.. include:: ../../../input.create_templates.dir/aquila_kato.tvel
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
