Template Kurtosis-based Waveform Check
**************************************

Kurtosis statistics is used to evaluate the simmetry of the time series neglecting from the pool of trimmed templates
waveforms that show high simmetry in the S-wave selected time window. It is supposed that low signal to noise ratios 
have low values of Kurtosis index. This contributes to exclude the template/channel from the estimation of the cross-correlation.
The routine avoids also signals ior glitches that are located at the beginning and at the end of the signal (please check carefully the code
before using it). 
Scipy Kurtosis is used as a standard routine (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kurtosis.html)

An example using the input_template_check is provided also showing accpted and removed waveforms. 

