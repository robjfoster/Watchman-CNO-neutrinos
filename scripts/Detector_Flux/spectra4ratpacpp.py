import numpy as np
import matplotlib.pyplot as plt
import pdb
import os
from scipy import integrate

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_pp = 5.98e10


# read data for the pp chain neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
pp = np.genfromtxt(dirname + "/../../data/pp.dat", delimiter="  ") # sns.ias.edu

#print(pp[:,2])

Epp1 = pp[:, 0]
fpp1 = pp[:, 1]*fi_pp
Epp2 = pp[:, 2]
fpp2 = pp[:, 3]*fi_pp
Epp3 = pp[:, 4]
fpp3 = pp[:, 5]*fi_pp
Epp4 = pp[:, 6]
fpp4 = pp[:, 7]*fi_pp

#print(Epp1,Epp2,Epp3,Epp4)

Epp = np.concatenate((Epp1,Epp2,Epp3,Epp4))
fpp = np.concatenate((fpp1,fpp2,fpp3,fpp4))

FreeElectrons = 3.34e32
nktons=1.32 #Approximate fiducial volume

me_si = 9.11e-31 # Electron mass kg
me=0.511


Epp_j = Epp*1.602e-13 #MeV to J

Gf_GeV = 1.1663788e-5 #GeV^-2
Gf_j = Gf_GeV/10**18/(1.602e-19)**2 # J^-2
Gf = Gf_GeV*1e-6 #MeV^-2
me_si = 9.11e-31 # Electron mass kg
me=0.511
hbar = 1.0545718e-34 # m^2 kg s^-1

sigma_0 = 8.806e-45 #cm**2 https://scholarworks.umass.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1109&context=dissertations_2

#sigma=sigma_0*Epp/me
#sigma_int=integrate.trapz(sigma,Epp)
#print(sigma_int)

sigma = 9.47e-45 * (Epp) #cm**2
#print(9.47e-45 * max(Enbinned))
sigma_int = integrate.trapz(sigma,Epp)
print(sigma_int)

Fbinned_int = integrate.trapz(fpp,Epp)
print(Fbinned_int)
rate = nktons*FreeElectrons*sigma_int*Fbinned_int
print("Rate = %.4e /s" % rate)
print("Rate = %.4e /day" % (rate*3600*24))

f = 6e11 #table 1.3 of Direct measurements of the pp solar neutrinos in Borexino

rate = FreeElectrons*sigma*6e11#Fbinned_int
rate_int = integrate.trapz(rate,Epp)
print(rate_int)

#with open("SolarNeutrinoRates.txt",'a+') as outfile:
#	outfile.write("Rate for pp in 1.32 kTon fiducial volume = %e per second\n" % rate)

#n_bins = len(Epp)

#write the combined spectrum to a file in the required format for ratpac
#with open(dirname + "/../../data/Detector_Flux/pp.ratdb",'w') as outfile:
#    outfile.write("{\n")
#    outfile.write('name: "SPECTRUM",\n')
#    outfile.write('index: "pp",\n')
#    outfile.write("valid_begin: [0,0],\n")
#    outfile.write("valid_End:[0,0],\n")
#    outfile.write("emin: %f\n" % Epp[0])
#    outfile.write("emax: %f\n" % Epp[-1])
#    outfile.write("spec_e: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % Epp[i])
#        else: outfile.write("%f],\n" % Epp[i])
#    outfile.write("spec_mag: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % Fpp[i])
#        else: outfile.write("%f],\n" % Fpp[i])
#    outfile.write("}")

