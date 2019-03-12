.. image:: pympa_logo1.png
    :width: 100px
    :align: left
    :alt: pympa_logo1.png
    :target: https://github.com/avuan/PyMPA37/releases



Main Program
============

PyMPA contains an efficient code for the detection of microseismicity starting from well located templates.
The code is stored on |github|, and is free to be cloned on your platform. It supports Python 2.7, 3.4, 3.5, 3.6, 3.7
releases and uses |ObsPy_link| for reading and writing seismic data, and for handling most
of the event metadata. Matched-filter correlations are calculated using ObsPy v. 1.2.0 |correlate_template| released on March
2019. 
Important: we recommend to use an updated version of ObsPy.
 
.. |github| raw:: html

    <a href="https://github.com/avuan/PyMPA37" target="_blank">github</a>

.. |ObsPy_link| raw:: html

  <a href="https://github.com/obspy/obspy/wiki" target="_blank">ObsPy</a>

.. |correlate_template| raw:: html

  <a href="https://docs.obspy.org/master/packages/autogen/obspy.signal.cross_correlation.correlate_template.html" target="_blank">correlate_template</a>


Main packages contains:

* :doc:`Template matching by using daily estimation of MAD and all the available channels <./sub/main.pympa>`;
* :doc:`Template matching by using daily estimation of MAD and a limited number of channels <./sub/main.pympa_channel_limit>`;
* :doc:`Template matching by using daily chunks (MAD estimated along the chunk duration) and a limited number of channels <./sub/main.pympa_chunks_channel_limit>`;

This package is written by the PyMPA developers, and is distributed under the LGPL GNU Licence, Copyright PyMPA developers 2019.


Contents:
---------
.. toctree::
   :maxdepth: 2

   main.pympa <./sub/main.pympa>
   main.pympa_channel_limit <./sub/main.pympa_channel_limit>
   main.pympa_chunks_channel_limit <./sub/main.pympa_chunks_channel_limit>
