Calculate Travel Times
**********************

Theoretical travel-time arrivals are calculated using the ObsPy port of the Java TauP Toolkit routines, see https://docs.obspy.org/packages/obspy.taup.html (Crotwell et al., 1999).
For using your own earth model see https://docs.obspy.org/packages/autogen/obspy.taup.taup_create.build_taup_model.html#obspy.taup.taup_create.build_taup_model
Model initialization is an expensive operation. Thus, make sure to do it only if necessary. 
ObsPy include custom built models can be initialized by specifying an absolute path to a model in ObsPy’s .npz model format instead of just a model name. See below for information on how to build a .npz model file.
Building an ObsPy model file from a “tvel” or “nd” file is easy.

An example of tvel model to be compiled by build_taup_model Obspy function:

.. include:: ../../../input.create_templates.dir/aquila_kato.tvel
   :start-line: 0
   :end-line: 10
   :literal:

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

References

Crotwell, H. P., T. J. Owens, and J. Ritsema (1999). The TauP Toolkit:
Flexible seismic travel-time and ray-path utilities, Seismol. Res. Lett.
70, 154–160.
