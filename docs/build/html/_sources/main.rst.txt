Main Program
============

PyMPA procedure provides the detection of microseismicity starting from well located templates by using a cross-correlation function from a network.
The code is stored on |github|, and can be freely cloned on your platform. It supports Python 2.7, 3.4, 3.5, 3.6, 3.7
releases and uses |ObsPy_link| for reading and writing seismic data, and for handling most
of the event metadata. Matched-filter correlations are calculated using a python normalised cross-correlation function or the
ObsPy v. 1.2.0 |correlate_template| released on April
2019. Detections can be also obtained using a single station and three channels by modifying the input parameters (thresholds etc..). 
This version is running on a single core and does not include multiprocessing routines. However, in the need of massive calculations for years and thousands
of templates, it could be easily implemented a script using SLURM or other schedulers to submit many jobs to the available processors.
   
Important: we recommend to use an updated version of ObsPy.
 
.. |github| raw:: html

    <a href="https://github.com/avuan/PyMPA37" target="_blank">github</a>

.. |ObsPy_link| raw:: html

  <a href="https://github.com/obspy/obspy/wiki" target="_blank">ObsPy</a>

.. |correlate_template| raw:: html

  <a href="https://docs.obspy.org/master/packages/autogen/obspy.signal.cross_correlation.correlate_template.html" target="_blank">correlate_template</a>


Main packages contains:

* :doc:`Template matching by using daily estimation of MAD and all the available channels <./sub/main.pympa>`;

This package is written by the PyMPA developers, and is distributed under the LGPL GNU Licence, Copyright PyMPA developers 2019.


Contents:
---------
.. toctree::
   :maxdepth: 2

   main.pympa <./sub/main.pympa>


.. image:: pympa_logo1.png
    :width: 100px
    :align: center
    :alt: pympa_logo1.png
    :target: https://github.com/avuan/PyMPA37/releases
