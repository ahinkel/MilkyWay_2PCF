The plotting code plots the data and makes any last minute calculations.  For example, argsPlotHist.py will take RR, DD, and DR histograms and compute the LS Estimator and then plot it.

I will use [input file here] for the various wildcards. Most of these are self explanatory, but I will detail one in particular:


TO RUN:
Python argsPlotHist.py [file containing RR DD and DR as columns, space delimited] [max separation distance to probe, assuming min is 0, in kpc.] [number of bins] [pdf file name to store plot figure as]


INPUT FILES:
File containing RR DD and DR as columns, space delimited (see 2PCF readme, output)


OTHER INPUTS:
User can specify the scale of the axes, axis titles, figure titles, etc. within the script.


OUTPUT:
The main output is of course the plot, specified above as [pdf file name to store plot figure as].  The script currently also prints out a couple of numbers to allow one to check if RR and DD have similar numbers of stars, which nominally they should.  Numbers close to one indicate similar histogram areas.  This part isn’t critical, but might be of interest, otherwise it can be commented out.  
