import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("AAS_prodrun/nofz_plots/nofz_hist_R7680_z1030NandS_phi180186.txt") #south then north

south = data[:,0]
north = data[:,1]

sNorm = 0.0 + np.sum(south)
nNorm = 0.0 + np.sum(north)

south = south/sNorm
north = north/nNorm

#reverse south data:
backwardSouth = np.zeros(len(south))
for i in range(0,len(south)):
    backwardSouth[i] = south[len(south) - i - 1]

#UPDATE THESE!!!!!!!!!!!!!!!!
zmin = -3.0 #kpc
zmax = 3.0 #kpc
Npoints = 200
z = np.linspace(zmin, zmax, num=2*Npoints)

#print sNorm, nNorm

#print south/sNorm

#paste together N and S:
NandS = np.concatenate((backwardSouth,north), axis=None)
#print NandS

#for fit, want to get rid of empty middle bins (|z|<0.2)
zForFit = np.array([z[0]])
NSforFit = np.array([NandS[0]])
for j in range(1,len(NandS)):
    if(NandS[j] != 0):
        zForFit = np.append(zForFit, z[j])
        NSforFit = np.append(NSforFit, NandS[j])

#########Chebyshev fit:
#x = z/zmax
deg = 64
#polynomialSouth = np.polynomial.Chebyshev.fit(z, south/sNorm, deg, domain=[0.2,2])
#polynomialNorth = np.polynomial.Chebyshev.fit(z, north/nNorm, deg, domain=[0.2,2])
#polynomial = np.polynomial.Chebyshev.fit(z, NandS, deg, domain=[-2,2])
polynomial = np.polynomial.Chebyshev.fit(zForFit, NSforFit, deg, domain=[zmin,zmax])
#UPDATE ABOVE DOMAIN!!!!!!!
#######################

#print polynomial 
coef = polynomial.coef
for i in range(deg):
    print coef[i]


#Quality of Fit:
"""
chiSq = 0.0;
for k in range(0,len(zForFit)):
    dif = np.abs(polynomial(zForFit[k]) - NSforFit[k])
    sig = np.sqrt(NSforFit[k])
    chiSq += (dif**2)/sig**2

dof = len(NSforFit) - deg
reducedChiSq = chiSq/dof
print "reduced Chi Sq: "
print reducedChiSq
"""

"""
moreZ = np.linspace(zmin,zmax,0.01)
plt.plot(zForFit, polynomial(zForFit),'k.')
plt.plot(moreZ, polynomial(moreZ),'r-')
plt.xlim(1,zmax)
plt.show()
"""

#plt.plot(z, south/sNorm, 'r.')
#plt.plot(z, north/nNorm, 'k.')
#plt.plot(z, NandS, 'k.')
plt.plot(zForFit, NSforFit, 'k.')
#plt.plot(z, polynomialSouth(z), 'm-')
plt.plot(zForFit, polynomial(zForFit), 'b-')

#plt.axhline(y=0,linestyle='-',color='k') 

plt.xlabel("z (kpc)")
plt.ylabel("normalized n(z) and Chebyshev fit ")
plt.xlim(zmin,zmax) #UPDATE SCALES!
plt.ylim(0, 0.035)
plt.title("normalized star counts fit (degree 64)")
plt.savefig("AAS_prodrun/nofz_plots/nofz_chebyshevfit_R7680_z1030_phi180186_deg64.pdf")
plt.show()

