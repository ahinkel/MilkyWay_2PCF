The code in this folder is used to compute the RR, DD, and DR histograms, and each script can compute this for a different coordinate separation.  This is the most time intensive portion of the analysis. Run time scales as N^2

I will use [input file here] for the various wildcards.


TO RUN:
g++ -O3 [one of the 2pcf scripts]  
./a.out  [input gaia data or mock galaxy data for RR histogram] [input gaia data or mock galaxy data for DD histogram] > [output file to output RR, DD, and DR to]


FILE INPUT:
For input, the galaxy files must be in Cylindrical coordinates, with columns ordered as: R, phi, z,   and the columns are space delimited.  As an exception, the  args3d_qOnly_xyzInput_2PCF.cpp code takes input files in cartesian coordinates: x, y, z,  also space delimited.


OTHER INPUTS:
In the code, the separation scales to probe are customizable (binMin and binMax, in kpc), as is the number of bins (numBins).  For example if one wanted to probe separation distances of 0pc to 750pc, with 100 bins, one would set: binMin = 0, binMax = 0.75, and numBins = 100.


OUTPUT:
A '.txt' file is generated and is a space delimited file with columns RR, DD, DR.



