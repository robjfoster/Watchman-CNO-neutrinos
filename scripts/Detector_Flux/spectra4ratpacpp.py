import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

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

Rfid = (10026.35-3080)/10 #Radius fiducial volume in cm
hfid = 2*Rfid #Height of fiducial volume in cm

theta = 1.3103 * np.pi/180 #Average Solar elevation at Boubly

Aeff = (np.cos(theta) * hfid * 2 * Rfid) + (np.sin(theta) * np.pi * Rfid**2) #Effective area of fiducial volume exposed to Sun

#fpp = np.pi * Rfid**2 * fpp #Flux in detector s^-1 MeV^-1

fpp = Aeff * fpp #Flux in detector s^-1 MeV^-1

FreeProtons = 0.668559 # 0.668559 * 10^32 Free protons per kton of water
nktons = 2 * np.pi * (hfid/100) * (Rfid/100)**2 /1000 # Number of ktons of water in detector. 1000 m^3 in 1 kton
TNU = FreeProtons * nktons #Using s^-1 not year^-1 as this is what watchmakers uses.

Fpp = fpp * TNU * Epp #Flux in terms of TNU (using s^-1)

n_bins = len(Epp)

#write the combined spectrum to a file in the required format for ratpac
with open(dirname + "/../../data/Detector_Flux/pp.ratdb",'w') as outfile:
    outfile.write("{\n")
    outfile.write('name: "SPECTRUM",\n')
    outfile.write('index: "pp",\n')
    outfile.write("valid_begin: [0,0],\n")
    outfile.write("valid_End:[0,0],\n")
    outfile.write("emin: %f\n" % Epp[0])
    outfile.write("emax: %f\n" % Epp[-1])
    outfile.write("spec_e: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % Epp[i])
        else: outfile.write("%f],\n" % Epp[i])
    outfile.write("spec_mag: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % Fpp[i])
        else: outfile.write("%f],\n" % Fpp[i])
    outfile.write("}")

