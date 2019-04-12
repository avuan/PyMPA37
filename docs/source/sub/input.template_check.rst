Template Kurtosis-based Waveform Check
**************************************

The input data quality of continuous waveforms is verified by testing daily gaps and overlapping streams and ensuring that a certain percentage of the requested time span is available. To obtain accurate results, input preparation is critical, and template waveforms should be carefully checked before running PyMPA. Filter frequency band selection depends on the best signal-to-noise ratio of templates and is related to the structural model and epicentral distances involved.
|Kurtosis| statistics is used to evaluate the simmetry of the time series neglecting from the pool of trimmed templates
waveforms that show high simmetry in the S-wave selected time window. It is supposed that low signal to noise ratios 
have low values of |Kurtosis| index. This contributes to exclude the template/channel from the estimation of the cross-correlation.
The routine avoids also signals or glitches that are located at the beginning and at the end of the signal (please check carefully the code before using it). 
Scipy |Kurtosis| is used as a standard routine. 

(see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kurtosis.html)

An example using the input_template_check is provided also showing accpted and removed waveforms. 


.. |Kurtosis| raw:: html

    <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kurtosis.html" target="_blank">Kurtosis</a>


