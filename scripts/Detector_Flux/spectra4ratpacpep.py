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

fi_pep = 7.98e8

Epep = 1.44 # Solar Neutrinos Spectroscopy 1704.06331

#EPEP=[Epep]*240
fPEP=np.linspace(0,fi_pep,240) #Create 240 points for PEP spectrum between min and max flux at fixed energy

EPEP=np.linspace(Epep-0.5*(240/10000000),Epep+0.5*(240/10000000),240)

FreeElectrons = 3.34e32
nktons=1.32 #Approximate fiducial volume

me_si = 9.11e-31 # Electron mass kg
me=0.511


EPEP_j = EPEP*1.602e-13 #MeV to J

Gf_GeV = 1.1663788e-5 #GeV^-2
Gf_j = Gf_GeV/10**18/(1.602e-19)**2 # J^-2
Gf = Gf_GeV*1e-6 #MeV^-2
me_si = 9.11e-31 # Electron mass kg
me=0.511
hbar = 1.0545718e-34 # m^2 kg s^-1

sigma_0 = 8.806e-45 #cm**2 https://scholarworks.umass.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1109&context=dissertations_2

#sigma=sigma_0*EPEP/me
#sigma_int=integrate.trapz(sigma,EPEP)
#print(sigma_int)

sigma = 9.47e-45 * (EPEP) #cm**2
#print(9.47e-45 * max(Enbinned))
sigma_int = integrate.trapz(sigma,EPEP)
print(sigma_int)

Fbinned_int = integrate.trapz(fPEP,EPEP)

rate = nktons*FreeElectrons*sigma_int*Fbinned_int
print("Rate = %.4e /s" % rate)
print("Rate = %.4e /day" % (rate*3600*24))

n_bins = len(EPEP)

energy_hist = np.genfromtxt(dirname + "/particle_energy_values_pep.txt", delimiter=",")

print(len(energy_hist))

(n,bins,patches)=plt.hist(energy_hist,n_bins)
plt.xlabel('Energy (MeV)')
plt.ylabel('Events')
plt.show()

print(max(fPEP))
print(Fbinned_int)

rate_plot = 3600*24*25*365.25*FreeElectrons*sigma*Fbinned_int#fmean#fbinned
rate_plot_unit=rate_plot*n[0]/(0.02*EPEP)/len(energy_hist)#/n_bins

fitted=gaussian_filter1d(rate_plot_unit, sigma=3)

plt.plot(EPEP,rate_plot_unit)
plt.plot(EPEP,fitted)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(0.6,7)
plt.show()

plt.plot(EPEP,fPEP)
plt.xscale('log')
plt.yscale('log')
plt.show()

#with open("SolarNeutrinoRates.txt",'a+') as outfile:
#	outfile.write("Rate for pep in 1.32 kTon fiducial volume = %e per second\n" % rate) #Actually roughly CNO rate

#n_bins = len(EPEP)

#write the combined spectrum to a file in the required format for ratpac
#with open(dirname + "/../../data/Detector_Flux/pep.ratdb",'w') as outfile:
#    outfile.write("{\n")
#    outfile.write('name: "SPECTRUM",\n')
#    outfile.write('index: "pep",\n')
#    outfile.write("valid_begin: [0,0],\n")
#    outfile.write("valid_end:[0,0],\n")
#    outfile.write("emin: %f\n" % Epep)
#    outfile.write("emax: %f\n" % Epep)
#    outfile.write("spec_e: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % EPEP[i])
#        else: outfile.write("%f],\n" % EPEP[i])
#    outfile.write("spec_mag: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % FPEP[i])
#        else: outfile.write("%f],\n" % FPEP[i])
#    outfile.write("}")

