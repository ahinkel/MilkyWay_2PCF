These scripts help to preprocess the data in various ways.  Most importantly, the “prepareCorrelationData” script set will apply cuts on the data via root, and output the smaller subset as a .txt file, space delimited, with columns R, phi, z.

Note makefile has warnings for unused vars -- can ignore.

TO RUN:
make -f Makefile.prepareCorrelationData
./prepareCorrelationData > [myFile2write2.txt]


FILE INPUTS:
Makefile, header file, and .C macro all need to be in same directory, and the root file is included on this google drive.


OTHER INPUTS:
Cuts are specified on line 67 of the code, to be implemented in root. Note all variable names in root are listed on line 41 of prepareCorrelationData.C 

These cuts will carve out the desired portion of the larger Gaia dataset.


OUTPUT:
A slice or wedge of Gaia data, in txt format  space delimited, with columns R, phi, z. Stored in [myFile2write2.txt]



