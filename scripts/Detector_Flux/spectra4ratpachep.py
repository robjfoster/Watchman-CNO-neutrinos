import numpy as np
import matplotlib.pyplot as plt
import pdb
import os
from scipy import integrate
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from scipy.ndimage.filters import gaussian_filter1d

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_hep = 5.28e3


# read data for the hep chain neutrinos. Flux in cm^-2 s^-1 MeV^-1 at Earth's surface, E in MeV
hep = np.genfromtxt(dirname + "/../../data/hep.dat", delimiter="  ") # sns.ias.edu

Ehep = hep[:, 0]
fhep = hep[:, 1]*fi_hep

FreeElectrons = 3.34e32
nktons=1.32 #Approximate fiducial volume

me_si = 9.11e-31 # Electron mass kg
me=0.511


Ehep_j = Ehep*1.602e-13 #MeV to J

Gf_GeV = 1.1663788e-5 #GeV^-2
Gf_j = Gf_GeV/10**18/(1.602e-19)**2 # J^-2
Gf = Gf_GeV*1e-6 #MeV^-2
me_si = 9.11e-31 # Electron mass kg
me=0.511
hbar = 1.0545718e-34 # m^2 kg s^-1

sigma_0 = 8.806e-45 #cm**2 https://scholarworks.umass.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1109&context=dissertations_2

#sigma=sigma_0*Ehep/me
#sigma_int=integrate.trapz(sigma,Ehep)
#print(sigma_int)

sigma = 9.47e-45 * (Ehep) #cm**2
#print(9.47e-45 * max(Enbinned))
sigma_int = integrate.trapz(sigma,Ehep)
print(sigma_int)

Fbinned_int = integrate.trapz(fhep,Ehep)

rate = nktons*FreeElectrons*sigma_int*Fbinned_int
print("Rate = %.4e /s" % rate)
print("Rate = %.4e /day" % (rate*3600*24))

n_bins = len(Ehep)

energy_hist = np.genfromtxt(dirname + "/particle_energy_values_hep.txt", delimiter=",")

print(len(energy_hist))

(n,bins,patches)=plt.hist(energy_hist,n_bins)
plt.xlabel('Energy (MeV)')
plt.ylabel('Events')
plt.show()

print(max(fhep))
print(Fbinned_int)

rate_plot = 3600*24*25*365.25*FreeElectrons*sigma*Fbinned_int#fmean#fbinned
rate_plot_unit=rate_plot*n[0]/(0.02*Ehep)/len(energy_hist)#/n_bins

fitted=gaussian_filter1d(rate_plot_unit, sigma=2)

plt.plot(Ehep,rate_plot_unit)
plt.plot(Ehep,fitted)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(0.6,7)
plt.show()

plt.plot(Ehep,fhep)
plt.xscale('log')
plt.yscale('log')
plt.show()


#with open("SolarNeutrinoRates.txt",'a+') as outfile:
#	outfile.write("Rate for hep in 1.32 kTon fiducial volume = %e per second\n" % rate)

#n_bins = len(Ehep)

#with open(dirname + "/../../data/Detector_Flux/hep.ratdb",'w') as outfile:
#    outfile.write("{\n")
#    outfile.write('name: "SPECTRUM",\n')
#    outfile.write('index: "hep",\n')
#    outfile.write("valid_begin: [0,0],\n")
#    outfile.write("valid_end:[0,0],\n")
#    outfile.write("emin: %f\n" % Ehep[0])
#    outfile.write("emax: %f\n" % Ehep[-1])
#    outfile.write("spec_e: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % Ehep[i])
#        else: outfile.write("%f],\n" % Ehep[i])
#    outfile.write("spec_mag: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % Fhep[i])
#        else: outfile.write("%f],\n" % Fhep[i])
#    outfile.write("}")

