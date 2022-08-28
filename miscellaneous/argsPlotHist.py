import sys
import numpy as np
import matplotlib.pyplot as plt

#TO RUN:
#python argsPlotHist.py <fileInputName> <maxDistance2Probe2> <numBins> <outputFileName>

inputFile = str(sys.argv[1])
maxR = float(sys.argv[2]); #kpc
numBins = int(sys.argv[3]);
fileForFigure = str(sys.argv[4]);

data = np.loadtxt(inputFile)

RR = data[:,0]
DD = data[:,1]
DR = data[:,2]

#delmethis part see if RR<-->DD swap matters (it shouldn't)
#temp = RR
#RR = DD
#DD = temp
#end delme

rrm1o2 = np.sum(RR) #number of entries in RR is r(r-1)/2
nnm1o2 = np.sum(DD) #number of entries in RR is n(n-1)/2
rn = np.sum(DR)     #number of entries in DR is r*n

RRerr = np.sqrt(RR)
DDerr = np.sqrt(DD)
DRerr = np.sqrt(DR)

d = np.linspace(0,maxR,num=numBins)

#careful with integer division / how do bignums work here?
A = rrm1o2/nnm1o2
B = 2.0*rrm1o2/rn

print A, B

#Correlators:
NaturalCorrelator = A * (DD.astype(float) / RR) - 1.0
NaturalCorrErr = (A/RR) * np.sqrt(DDerr**2 + RRerr**2 * (DD.astype(float)/RR)**2)
#Error calcs updated and are correct, 2020-Nov-19

LSCorr = 1.0 + A*(DD.astype(float))/RR - B*(DR.astype(float))/RR
LSCorrErr = (1.0/RR)*np.sqrt(A*A*DDerr**2 + B*B*DRerr**2 + RRerr**2*((A*DD-B*DR)/RR.astype(float))**2)


#plt.errorbar(d, NaturalCorrelator, NaturalCorrErr, capsize=1, marker='.')
#plt.plot(d, AltNatCorr, 'g--')
plt.errorbar(d, LSCorr, LSCorrErr, capsize=1, marker='.')
plt.axhline(y=0,linestyle='-',color='k') 
plt.xlabel("Separation Distance, x-only (kpc)")
#plt.ylabel("2PCF ")
plt.ylabel(r'$\xi$')
#plt.xlim(0,0.25)
plt.xlim(0,maxR) #UPDATE SCALES!
#plt.xlim(0,0.5)
plt.ylim(-0.005, 0.005)
#plt.yscale('log')
#plt.xscale('log')
#plt.title("Model vs. Model Example")
#plt.title("Halo Study: Model vs. Gaia")
#plt.title("Gaia (North) vs. Gaia (South)")
plt.tight_layout() #avoids the y-axis being cut off
plt.savefig(fileForFigure)
plt.show()
