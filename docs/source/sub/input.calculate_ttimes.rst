Calculate Travel Times
**********************

Needed files:

- Events in a catalog: e.g. templates.zmap (quakeml or zmap format) see ObsPy for format

- Suitable velocity model for computing travel times

- Station inventory (format consistent with ObsPy read_inventory routine see https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.read_inventory.html)

- Days to process: one column file including days to process e.g. lista1

- Set parameters: e.g. times.par

.. include:: ../../../input.calculate_ttimes.dir/times.par
   :literal:

- Input directory i.e. ./template where trimmed templates are found

- Output dir i.e. ./ttimes (find moveout times from different channels used to synchronize cross-correlation functions)
