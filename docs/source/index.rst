.. image:: pympa_logo.png
    :width: 600px
    :align: left
    :alt: pympa_logo.png
    :target: https://github.com/avuan/PyMPA37/releases

PyMPA 

Python Matching Phase Algorithm
==========

A Python package for the detection of seismicity based on templates.
PyMPA contains an efficient code for the detection of microseismicity starting from well located templates.
The software package PyMPA is an open source seismological software. It consists of some separate utilities for data preparation, the main program, and post-processing tools to obtain a catalog and verify events. PyMPA is designed to detect microseismicity from the cross-correlation of continuous data and templates.


:doc:`input_preparation </submodules/core.input_preparation>` :doc:`detection </submodules/core.detection>`
:doc:`postprocessing </submodules/core.postprocessing>`, and verification :doc:`verify </submodules/core.verify>`.  

The code is stored on |github|, and is free to be cloned on your platform. It supports Python 2.7, 3.4, 3.5, 3.6, 3.7releases and uses |ObsPy_link| for reading and writing seismic data, and for handling most
of the event metadata. Matched-filter correlations are calculated using |correlate_template|.
 
.. |github| raw:: html

    <a href="https://github.com/avuan/PyMPA37" target="_blank">github</a>

.. |ObsPy_link| raw:: html

  <a href="https://github.com/obspy/obspy/wiki" target="_blank">ObsPy</a>

.. |correlate_template| raw:: html

  <a href="https://docs.obspy.org/master/packages/autogen/obspy.signal.cross_correlation.correlate_template.html" target="_blank">correlate_template</a>

This package contains:

* :doc:`Routines for seismic data preparation </submodules/utils.preparation>`;
* :doc:`Kurtosis based template verification </submodules/utils.kurtosis>`;
* :doc:`Postprocessing routines </submodules/utils.postprocessing>`;
* :doc:`Detection visual verification </submodules/utils.verification>`;

This package is written by the PyMPA developers, and is 
distributed under the LGPL GNU Licence, Copyright PyMPA
developers 2019.

Citation
--------
If you use this package in your work, please cite the following papers:

Vuan A., Sugan M., Amati G., Kato A., 2017 - Improving the Detection of Low-Magnitude Seismicity Preceding the Mw 6.3 L'Aquila Earthquake: Development of a Scalable Code Based on the Cross-Correlation of Template Earthquakes, BSSA https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/525813/improving-the-detection-of-low-magnitude?redirectedFrom=fulltext

Vuan A., Sugan M., Chiaraluce L., Di Stefano R., 2017 - Loading rate variations along a mid-crustal shear zone preceding the MW6.0 earthquake of the 24th of August 2016 in Central Italy, Geophysical Research Letters http://onlinelibrary.wiley.com/doi/10.1002/2017GL076223/full

Sugan, M., Vuan, A., Kato, A., Massa, M., & Amati, G. (2019). Seismic evidence of an early afterslip during the 2012 sequence in Emilia (Italy). Geophysical Research Letters, 46, 625–635. https://doi.org/10.1029/2018GL079617

Contents:
---------

.. toctree::
   :numbered:
   :maxdepth: 2

   intro
   installation
   updates
   tutorial
   core
   utils
