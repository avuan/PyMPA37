.. image:: pympa_logo1.png
    :width: 100px
    :align: left
    :alt: pympa_logo1.png
    :target: https://github.com/avuan/PyMPA37/docs

Python Matching Phase Algorithm 
===============================

A Python package for the detection of seismicity based on templates.
PyMPA contains an efficient code for the detection of microseismicity starting from well located templates.
The software package PyMPA is an open source seismological software. It consists of some separate utilities for input preparation, the main program, and output post-processing tools to obtain a catalog and verify events.  PyMPA is designed to detect microseismicity from the cross-correlation of continuous data and templates.

The code is stored on |github|, and is free to be cloned on your platform. It supports Python 2.7, 3.4, 3.5, 3.6, 3.7 releases and uses |ObsPy_link| for reading and writing seismic data, and for handling most
of the event metadata. Matched-filter correlations are calculated using ObsPy v. 1.2.0 |correlate_template| released on March
2019. 
Important: we recommend to use an updated version of ObsPy.
 
.. |github| raw:: html

    <a href="https://github.com/avuan/PyMPA37" target="_blank">github</a>

.. |ObsPy_link| raw:: html

  <a href="https://github.com/obspy/obspy/wiki" target="_blank">ObsPy</a>

.. |correlate_template| raw:: html

  <a href="https://docs.obspy.org/master/packages/autogen/obspy.signal.cross_correlation.correlate_template.html" target="_blank">correlate_template</a>

.. image:: pympa_logo.png
    :width: 600px
    :align: center
    :alt: pympa_logo.png
    :target: https://github.com/avuan/PyMPA37/docs

This package contains:

- :doc:`Routines for downloading data from eida servers <./sub/input.download_data>`;
- :doc:`Routines for creating and trimming templates <./sub/input.create_templates>`;
- :doc:`Routines for calculating moveout time for synchronization <./sub/input.calculate_ttimes>`;
- :doc:`Kurtosis based template verification <./sub/input.template_check>`;
- :doc:`Template matching by using daily estimation of MAD and all the available channels <./sub/main.pympa>`;
- :doc:`Postprocessing routines <./sub/output.process_detections>`;
- :doc:`Visual verification of detections <./sub/output.verify_detection>`;

This package is written by the PyMPA developers, and is distributed under the LGPL GNU Licence, Copyright PyMPA developers 2019.


Acknowledgements
----------------
The software development was partially funded by a joint research project within the
executive program of scientific and technological cooperation between Italy
and Japan for the period 2013–2015. Additional funds for software development
come from the project “Seismology and Earthquake Engineering
Research Infrastructure Alliance for Europe” (SERA), responding to the priorities
identified in the call INFRAIA-01-2016-2017 Research Infrastructure
for Earthquake Hazard. We thank Monica Sugan for the extensive testing of the codes and Aitaro Kato 
at the Earthquake Research Institute (ERI) in Tokyo for fruitful
discussions.
The authors also wish to thank the ObsPy community for the continuous
support and constant development of related libraries.


Citation
--------
If you use this package in your work, please cite the following papers:

Vuan A., Sugan M., Amati G., Kato A., 2017 - Improving the Detection of Low-Magnitude Seismicity Preceding the Mw 6.3 L'Aquila Earthquake: Development of a Scalable Code Based on the Cross-Correlation of Template Earthquakes, BSSA https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/525813/improving-the-detection-of-low-magnitude?redirectedFrom=fulltext

Vuan A., Sugan M., Chiaraluce L., Di Stefano R., 2017 - Loading rate variations along a mid-crustal shear zone preceding the MW6.0 earthquake of the 24th of August 2016 in Central Italy, Geophysical Research Letters http://onlinelibrary.wiley.com/doi/10.1002/2017GL076223/full

Sugan, M., Vuan, A., Kato, A., Massa, M., & Amati, G. (2019). Seismic evidence of an early afterslip during the 2012 sequence in Emilia (Italy). Geophysical Research Letters, 46, 625–635. https://doi.org/10.1029/2018GL079617

Contents:
---------
.. toctree::
   :maxdepth: 2 
   :numbered:

   intro
   installation
   tutorial
   input
   main
   output
