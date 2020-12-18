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

fi_be2 = 1.44e9

# read data for the second be-7 set of neutrinos. Flux in cm^-2 s^-1 MeV^-1 at Earth's surface, E in MeV
be2 = np.genfromtxt(dirname + "/../../data/be7_2.dat", delimiter="  ") # sns.ias.edu
Ebe2 = (be2[:, 0]/1000) + .8613 # KeV to MeV and addition to return to energy
fbe2 = be2[:, 1]*fi_be2

FreeElectrons = 3.34e32
nktons=1.32 #Approximate fiducial volume

me_si = 9.11e-31 # Electron mass kg
me=0.511


Ebe2_j = Ebe2*1.602e-13 #MeV to J

Gf_GeV = 1.1663788e-5 #GeV^-2
Gf_j = Gf_GeV/10**18/(1.602e-19)**2 # J^-2
Gf = Gf_GeV*1e-6 #MeV^-2
me_si = 9.11e-31 # Electron mass kg
me=0.511
hbar = 1.0545718e-34 # m^2 kg s^-1

sigma_0 = 8.806e-45 #cm**2 https://scholarworks.umass.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1109&context=dissertations_2

#sigma=sigma_0*Ebe2/me
#sigma_int=integrate.trapz(sigma,Ebe2)
#print(sigma_int)

sigma = 9.47e-45 * (Ebe2) #cm**2
#print(9.47e-45 * max(Enbinned))
sigma_int = integrate.trapz(sigma,Ebe2)
print(sigma_int)

Fbinned_int = integrate.trapz(fbe2,Ebe2)

rate = nktons*FreeElectrons*sigma_int*Fbinned_int
print("Rate = %.4e /s" % rate)
print("Rate = %.4e /day" % (rate*3600*24))

n_bins = len(Ebe2)

energy_hist = np.genfromtxt(dirname + "/particle_energy_values_7Be_2.txt", delimiter=",")

print(len(energy_hist))

(n,bins,patches)=plt.hist(energy_hist,n_bins)
plt.xlabel('Energy (MeV)')
plt.ylabel('Events')
plt.show()

print(max(fbe2))
print(Fbinned_int)

rate_plot = 3600*24*25*365.25*FreeElectrons*sigma*Fbinned_int#fmean#fbinned
rate_plot_unit=rate_plot*n[0]/(0.02*Ebe2)/len(energy_hist)#/n_bins

fitted=gaussian_filter1d(rate_plot_unit, sigma=2)

plt.plot(Ebe2,rate_plot_unit)
plt.plot(Ebe2,fitted)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(0.6,7)
plt.show()

plt.plot(Ebe2,fbe2)
plt.xscale('log')
plt.yscale('log')
plt.show()

f = ((50+45.6)/2)*1e8 #max(fbe2) table 1.3 of Direct measurements of the pp solar neutrinos in Borexino

rate = f*FreeElectrons*sigma
rate_int = integrate.trapz(rate,Ebe2)
print(rate_int)


#with open("SolarNeutrinoRates.txt",'a+') as outfile:
#	outfile.write("Rate for 7Be_2 in 1.32 kTon fiducial volume = %e per second\n" % rate)

#n_bins = len(Ebe2)

#write the combined spectrum to a file in the required format for ratpac
#with open(dirname + "/../../data/Detector_Flux/be2.ratdb",'w') as outfile:
#    outfile.write("{\n")
#    outfile.write('name: "SPECTRUM",\n')
#    outfile.write('index: "7Be_2",\n')
#    outfile.write("valid_begin: [0,0],\n")
#    outfile.write("valid_end:[0,0],\n")
#    outfile.write("emin: %f\n" % Ebe2[0])
#    outfile.write("emax: %f\n" % Ebe2[-1])
#    outfile.write("spec_e: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % Ebe2[i])
#        else: outfile.write("%f],\n" % Ebe2[i])
#    outfile.write("spec_mag: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % Fbe2[i])
#        else: outfile.write("%f],\n" % Fbe2[i])
#    outfile.write("}")

