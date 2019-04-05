Process PyMPA Output
--------------------

From the defined template by cross-correlation for a selected day is obtained a list of possible detections.
Since time overlapping detections could be found also by different templates, the procedure allows a search for 
the template able to detect the events with the highest threshold value that is also related with the highest 
average cross-correlation value for the used network. From the main program, many cat files are daily sorted
and grouped for a scan of the events showing the best detections. A bash script is used to collect data and create the input for the procedure
and a python script is used to filter the result on the basis of the fiter.par parameters used. The filter.par defines
some additional filtering to detections to possibly overcome a visual validation of the new detections.  
