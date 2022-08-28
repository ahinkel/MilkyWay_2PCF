import numpy as np
import matplotlib.pyplot as plt

print "Don't forget to update length scales."
maxR = 0.40; #kpc
numBins = 200;

#data = np.loadtxt("AAS_prodrun/GAIA_GAIA_correlationData/zSun_zRefl/corrdata_zonly_R7476_z02151985NandSrefl_phi179180.txt")
#data = np.loadtxt("AAS_prodrun/GAIA_GAIA_correlationData/z_Refl/yonly_Corr/y350bin_0000pc_0350pc_R8084_z0207NvS_phi177180.txt")
data = np.loadtxt("erkalXYZ_corr_x200bin_0000pc_0400pc_R7882_z0220N_phi180181_seeds779v144_60kSteps.txt")
#data = np.loadtxt("AAS_prodrun/correlationData/zSun_zOnly/corrdata_R8082_z02151985N_phi179180_seed1993.txt")
#data = np.loadtxt("AAS_prodrun/GAIA_GAIA_correlationData/phi_Refl/R_only_Corr/Ronly_corr_R7090_z0203S_phi179180vs180181.txt")
#data = np.loadtxt("AAS_prodrun/GAIA_GAIA_correlationData/phi_Refl/phi_only/refl_corr_R8384_z0210N_phi174180_andrefl.txt")
#data = np.loadtxt("AAS_prodrun/correlationData/chebyshevZonly/800bins_0000pc_1000pc_R7684_z1230N_phi178182_seed442_deg32_lmcsmc.txt") #dont forget b/lmc/smc cuts
#data = np.loadtxt("AAS_prodrun/correlationData/phiOnly_chebyshev/corr400bins_0000pc_0400pc_R7684_z0508N_phi178182_seed4231_deg32_lmcsmc.txt")
#data = np.loadtxt("AAS_prodrun/correlationData/R_only_chebyshev/corr400bins_0000pc_0400pc_R7684_z0508N_phi178182_seed4231_deg32_lmcsmc.txt")
#data = np.loadtxt("AAS_prodrun/zSun_test/corrdata_R8486_z02150785_phi179181_zsun15pc_mockmock.txt")


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

#AltNatCorr = (DD.astype(float)/altRR) - 1.0
#print AltNatCorr
#print NaturalCorrErr
#print LSCorrErr

#plt.errorbar(d, NaturalCorrelator, NaturalCorrErr, capsize=1, marker='.')
#plt.plot(d, AltNatCorr, 'g--')
plt.errorbar(d, LSCorr, LSCorrErr, capsize=1, marker='.')
plt.axhline(y=0,linestyle='-',color='k') 
plt.xlabel("Separation Distance, x-only (kpc)")
plt.ylabel("2PCF ")
plt.xlim(0,maxR) #UPDATE SCALES!
plt.ylim(-0.005, 0.005)
#plt.yscale('log')
#plt.xscale('log')
plt.title("Model vs. Model")
#plt.savefig("AAS_prodrun/plots/phiOnly_chebyshev/400bins_0000pc_0400pc_R7684_z0508N_phi178182_seed4231_deg32_lmcsmc.pdf")
#plt.savefig("AAS_prodrun/GAIA_GAIA_plots/z_refl_plots/yonly_corr/y350bin_0000pc_0350pc_R8084_z0207NvS_phi177180.pdf")
plt.savefig("fig_erkalXYZ_corr_x200bin_0000pc_0400pc_R7882_z0220N_phi180181_seeds779v144_60kSteps.pdf")
plt.show()
