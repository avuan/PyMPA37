PyMPA
===

**a software package for phase match filtering**

![Screenshot](screenshot.png)

### Abstract

The software package PyMPA is an open source seismological software. It consists of some separate utilities for data preparation, the main program pympa38.mac.py, and post-processing tools to obtain a catalog and verify events. PyMPA is designed to detect microseismicity from the cross-correlation of continuous data and templates.

### Motivation and significance

The software was developed in the framework of a bilateral project between Italy and Japan 2013-2015 "Detection of slow slip events and space-time cluster analysis of low magnitude tremors preceding the 2009 L'Aquila earthquake and the 2012 Emilia seismic sequence in Italy", and was co-funded by Ministero Affari Esteri (Italy) - Scientific and Technological Collaboration Executive Program. Further improvements to the code are planned within Project H2020 - SERA 2017-2019 http://www.sera-eu.org/. 


### System requirements

Python installed to run the program (version 2.6 or more) is required [http://python.org]
Python 3 with Obspy installed (since version 0.10.1) is supported.
Installing Obspy implies Numpy, Scipy and Matplotlib libraries. For details on installing Obspy, please, consult the [official guide] (https://github.com/obspy/obspy/wiki).

### References
Vuan A., Sugan M., Amati G., Kato A., 2017 - Improving the Detection of Low-Magnitude Seismicity Preceding the Mw 6.3 L'Aquila Earthquake: Development of a Scalable Code Based on the Cross-Correlation of Template Earthquakes, BSSA
https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/525813/improving-the-detection-of-low-magnitude?redirectedFrom=fulltext

Vuan A., Sugan M., Chiaraluce L., Di Stefano R., 2017 - Loading rate variations along a mid-crustal shear zone preceding the MW6.0 earthquake of the 24th of August 2016 in Central Italy, Geophysical Research Letters http://onlinelibrary.wiley.com/doi/10.1002/2017GL076223/full
