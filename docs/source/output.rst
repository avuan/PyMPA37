.. image:: pympa_logo1.png
    :width: 100px
    :align: left
    :alt: pympa_logo1.png
    :target: https://github.com/avuan/PyMPA37/releases



Creating an output catalog and verify detections
================================================

    It consists of some separate utilities in |github|, post-processing tools to obtain a catalog and verify events.
    Many events could be correlated to more than one template in a narrow time window. A fixed time window length can be selected,
    and within each, the template for which the normalized cross-correlation coefficient is the greates provides the event location and data to determine the magnitude.
    This process is run by <./sub/output.process_detections> in to two steps,
    by using the last event origin time as a reference to set the next time window scrutinized. The final catalog should be verified by visual inspection
    for a number of sampled detections. Generally, we proceed by verification of events having low thresholds to understand a safe value to validate the catalog.
    The routine <./sub/output.verify_detection> creates graphs of time windows where continuous data and trimmed templates are plotted with info grasped from channel by channel cross-correlation process.


Important: we recommend to use an updated version of |ObsPy_Link|.
 
.. |github| raw:: html

  <a href="https://github.com/avuan/PyMPA37" target="_blank">github</a>

.. |ObsPy_link| raw:: html

  <a href="https://github.com/obspy/obspy/wiki" target="_blank">ObsPy</a>


These utilities contains:

* :doc:`Postprocessing routines <./sub/output.process_detections>`;
* :doc:`Visual verification of detections <./sub/output.verify_detection>`;


This package is written by the PyMPA developers, and is distributed under the LGPL GNU Licence, Copyright PyMPA developers 2019.

Contents:
---------
.. toctree::
   :maxdepth: 2

   output.process_detections <./sub/output.process_detections>
   output.verify_detection <./sub/output.verify_detection>

