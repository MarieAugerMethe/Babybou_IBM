Babybou_IBM
===========

Identifying parturition and neonate survival using movement data. 
This is the R code for the individual based method to identifying parturition and neonate survival from the movement data of female woodland caribou. 
The code was originally published in DeMars et al. (2013), and updated in Gharajehdaghipour et. al (2025). 
The code and original paper and supplementary materials are provided.

*Details on the update*:<br> 
In the updated version, break points representing parturition and neonate mortality events can occur at any time, including during gaps in the telemetry data. 
In contrast, in the previous version (see v1.0.0 and DeMars et al. 2013),
break points could only occur at time points when step length data are available. 
Therefore, NAs (indicating gaps in the telemetry data) are allowed in the numeric vector containing step lengths (*SL*), 
and the integer vector that identifies the time of each step length (*ti*).
*Important:* To be able to identify break points during gaps, you must create time series of *SL* and *ti* that are at regular time intervals.
Only the `IMB.R` function was modified, but similar changes could be made to `IMB_full_likelihood.R`.