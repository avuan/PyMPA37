.. image:: pympa_logo1.png
    :width: 100px
    :align: left
    :alt: pympa_logo1.png
    :target: https://github.com/avuan/PyMPA37/releases

PyMPA installation
=======================

PyMPA is a pure Python package. It runs after the installation of a virtual
environment with Numpy, Scipy, Matplotlib and Obspy libraries.  
Some C extensions of ObsPy toolkit are also used and Bottleneck libraries.
Bottleneck is a set of functions inspired by NumPy and SciPy, but written in 
Cython with high performance in mind. Bottleneck provides separate Cython 
functions for each combination of array dimensions, axis, and data type.

We heavily recommend installing ObsPy using conda because this will:

 * make your life easier;
 * separate your install from your system default Python, meaning you can
   experiment without breaking your OS;
 * ensure that modules are compiled correctly

If you do not have either a miniconda or anaconda installation you can follow
the |conda-install| instructions.

.. |conda-install| raw:: html

  <a href="https://www.anaconda.com/distribution/" target="_blank">conda-install</a>

If you do not already have a conda environment we recommend creating one
with the following:

.. code-block:: bash

    conda create -n obspy python=3.6
    conda activate obspy
    conda install obspy
    conda install bottleneck

For installing PyMPA you can simply clone the git repository and try it within your obspy env:

.. code-block:: bash
    
    git clone git@github.com:avuan/PyMPA37.git
