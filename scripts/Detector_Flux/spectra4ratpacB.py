import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_b = 2.78e6

# read data for the boron-8 neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
b8 = np.genfromtxt(dirname + "/../../data/b8.dat", delimiter="       ") # sns.ias.edu
Eb = b8[:, 0]
fb = b8[:, 1]*fi_b # b8.dat not delimited properly (replaced number of spaces in file to match)

Rfid = (10026.35-3080)/10 #Radius fiducial volume in cm
hfid = 2*Rfid #Height of fiducial volume in cm

theta = 1.3103 * np.pi/180 #Average Solar elevation at Boubly

Aeff = (np.cos(theta) * hfid * 2 * Rfid) + (np.sin(theta) * np.pi * Rfid**2) #Effective area of fiducial volume exposed to Sun

#fb = np.pi * Rfid**2 * b8[:, 1]*fi_b

fb = Aeff * b8[:, 1]*fi_b #s^-1 MeV^-1

FreeProtons = 0.668559 # 0.668559 * 10^32 Free protons per kton of water
nktons = 2 * np.pi * (hfid/100) * (Rfid/100)**2 /1000 # Number of ktons of water in detector. 1000 m^3 in 1 kton
TNU = FreeProtons * nktons #Using s^-1 not year^-1 as this is what watchmakers uses.

Fb = fb * TNU * Eb #Flux in terms of TNU (using s^-1)

n_bins = len(Eb)

#write the combined spectrum to a file in the required format for ratpac
with open(dirname + "/../../data/Detector_Flux/b8.ratdb",'w') as outfile:
    outfile.write("{\n")
    outfile.write('name: "SPECTRUM",\n')
    outfile.write('index: "B8",\n')
    outfile.write("valid_begin: [0,0],\n")
    outfile.write("valid_end:[0,0],\n")
    outfile.write("emin: %f\n" % Eb[0])
    outfile.write("emax: %f\n" % Eb[-1])
    outfile.write("spec_e: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % Eb[i])#, "%s" % "d,")
        else: outfile.write("%f],\n" % Eb[i])
    outfile.write("spec_mag: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % Fb[i])
        else: outfile.write("%f],\n" % Fb[i])
    outfile.write("}")

