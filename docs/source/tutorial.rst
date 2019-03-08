.. image:: pympa_logo1.png
    :width: 100px
    :align: left
    :alt: pympa_logo1.png
    :target: https://github.com/avuan/PyMPA37/releases

Tutorial 
======================================

This tutorial is designed to give you an overview of the capabilities and
implementation of the PyMPA Python package.

Downloading Seismological Data
------------------------------
To download seismological data from EIDA (European Integrated Data Archive) servers: 
Waveform broad band data from european broad band seismic stations are available 
from many others European Institutions. To download seismological data from EIDA (European Integrated Data Archive) servers
and inventory data in STATIONXML format many examples can be found in ObsPy examples.

PyMPA requires, continuous data and stations inventories.
EIDA servers can easily release continuous data from permanent networks and the corresponding 
inventories. The examples in the subdirectory input.download_data.dir show the python scripts
that allows the download. 

In the case your data come from other sources, PyMPA through ObsPy libraries
is able to manage most of the seismological data formats (MSEED, SAC, SEISAN, SEGY, etc..). 
An inventory data file including station information needs to be created by modifying an existing 
StationXML file and read that with ObsPy to perform the rest.

PyMPA does not use databases and prefers to store single channel daily continuous data in archieves.

(:doc:`download_data </submodules/input.download_data>`) download data

Create Templates
----------------
(:doc:`create_template </submodules/input.create_template>`) create templates

Needed files
1. Events in a catalog: e.g. templates.zmap (quakeml or zmap formats) see ObsPy for admitted formats
   cat templates.zmap
   
2. Days to process: one column file including days to process e.g. lista1
   cat lista1
   120103
   120104

3. Set parameters: e.g. trim.par
   cat trim.par
   #Line 1 -- list of stations
   #Line 2 -- list of channels
   #Line 3 -- list of networks
   #Line 4 -- Lowpass frequency
   #Line 5 -- Highpass frequency
   #Line 6 -- Trimmed Time before S-wave
   #Line 7 -- Trimmed Time after S-wave
   #Line 8 -- UTC precision
   #Line 9 -- Continuous data dir
   #Line 10 -- Template data dir
   #Line 11 -- Processing days list 
   #Line 12 -- Zmap catalog
   #Line 13 -- Starting template
   #Line 14 -- Stopping template
   #Line 15 -- Taup Model
   AQU CAMP CERT FAGN FIAM GUAR INTR MNS NRCA TERO
   BHE BHN BHZ
   IV MN
   2.0
   8.0
   2.5
   2.5
   6
   24h
   template
   lista1
   templates.zmap
   26
   27
   aquila_kato
   
Directories
-----------
./24h = include daily continuous time series 
./template = output trimmed templates

Velocity Model use to compute travel times
------------------------------------------
cat aquila_kato.tvel
aquila 
depth P vel. S vel. density older density
     0.000      3.7500      2.1650      2.4500 
     1.500      3.7500      2.1710      2.4500
     1.510      4.9400      2.8520      2.7800
     4.510      4.9400      2.8580      2.7800
     4.520      6.0100      3.2790      2.7600
    14.520      6.0100      3.2850      2.7600
    14.530      5.5500      3.3950      2.9100
    29.530      5.5500      3.4010      2.9100
    29.540      5.8800      4.0990      3.1000
    43.540      5.8800      4.1050      3.1000
    43.550      5.8800      4.5610      3.1000
    57.500      5.8800      3.3600      3.1000
    57.500      7.1100      4.0100      3.2300
    93.000      7.1100      4.0100      3.2300
    93.000      7.1000      3.9900      3.3000



Check Template Quality
----------------------
(:doc:`template_check </submodules/input.template_check>`) select good templates

Calculate Travel Times
----------------------
(:doc:`calculate_ttimes </submodules/input.calculate_ttimes>`) calculate travel times

Running PyMPA
-------------
(:doc:`main.pympa_chunks_channel_limit </submodules/main.pympa_chunks_channel_limit>`) run pympa

Output Processing
-----------------
(:doc:`output.process_detections </submodules/output.process_detections>`) controls multiple detections in short time windows

Verify Detections
-----------------
(:doc:`output.verify_detection </submodules/output.verify_detection>`) visual verification of events 




References
----------

Shelly, D. R., G. C. Beroza, and S. Ide (2007). Non-volcanic tremor and low
frequency earthquake swarms, Nature 446, 305–307.

Peng, Z., and P. Zhao (2009). Migration of early aftershocks following the
2004 Parkfield earthquake, Nature Geosci. 2, 877–881.

Yang, H., L. Zhu, and R. Chu (2009). Fault-plane determination of the
18 April 2008 Mount Carmel, Illinois, earthquake by detecting and
relocating aftershocks, Bull. Seismol. Soc. Am. 99, 3413–3420.

Kato, A., K. Obara, T. Igarashi, H. Tsuruoka, S. Nakagawa, and N. Hirata
(2012). Propagation of slow slip leading up to the 2011 Mw 9.0
Tohoku-Oki earthquake, Science 335, 705–708.

Zhang, M., and L. Wen (2015). An effective method for small event detection:
Match and locate (M&L), Geophys. J. Int. 200, 1523–1537.

Krischer, L., T. Megies, R. Barsch, M. Beyreuther, T. Lecocq, C. Caudron,
and J. Wassermann (2015). ObsPy: A bridge for seismology into the
scientific Python ecosystem, Comput. Sci. Discov. 8, no. 1, 014003,
doi: 10.1088/1749-4699/8/1/014003.
