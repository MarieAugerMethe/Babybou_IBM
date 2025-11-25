Babybou_IBM
===========

Identifying parturition and neonate survival using movement data. This is the R code for the individual based method to identifying parturition and neonate survival from the movement data of female woodland caribou. The code was originally published in DeMars et al. (2013), and updated (current version) in Gharajehdaghipour et. al (2025). The code and the papers and supplementary materials are provided.

*Details on the update*:<br> 
In the updated version presented here, break points representing parturition and neonate mortality events can occur at any time (i.e., including during gaps in the telemetry data), while the previous version published in Demars et al. (2013) only considered times where step length data are available. Therefore, NAs (indicating gaps in the telemetry data) are allowed in the numeric vector containing step lengths (*SL*), and the integer vector that identifies the time of each step length (*ti*).
