import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_be1 = 5e8

# read data for the first set of be-7 neutrinos. Flux in cm^-2 s^-1 MeV^-1 at Earth's surface, E in MeV
be1 = np.genfromtxt(dirname + "/../../data/be7_1.dat", delimiter="  ") # sns.ias.edu
Ebe1 = (be1[:, 0]/1000) + 0.3843 # KeV to MeV and addition to return to energy
fbe1 = be1[:, 1]*fi_be1

Rfid = (10026.35-3080)/10 #Radius fiducial volume in cm
hfid = 2*Rfid #Height of fiducial volume in cm

theta = 1.3103 * np.pi/180 #Average Solar elevation at Boubly

Aeff = (np.cos(theta) * hfid * 2 * Rfid) + (np.sin(theta) * np.pi * Rfid**2) #Effective area of fiducial volume exposed to Sun

#fbe1 = np.pi * Rfid**2 * be1[:, 1]*fi_be1 #Flux in detector s^-1 MeV^-1

fbe1 = Aeff * be1[:, 1]*fi_be1 #Flux in detector s^-1 MeV^-1

n_bins = len(Ebe1)

#write the combined spectrum to a file in the required format for ratpac
with open(dirname + "/../../data/Detector_Flux/be1.ratdb",'w') as outfile:
    outfile.write("{\n")
    outfile.write('name: "SPECTRUM",\n')
    outfile.write('index: "7Be_1",\n')
    outfile.write("valid_begin: [0,0],\n")
    outfile.write("valid_end:[0,0],\n")
    outfile.write("emin: %f\n" % Ebe1[0])
    outfile.write("emax: %f\n" % Ebe1[-1])
    outfile.write("spec_e: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % Ebe1[i])
        else: outfile.write("%f],\n" % Ebe1[i])
    outfile.write("spec_mag: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % fbe1[i])
        else: outfile.write("%f],\n" % fbe1[i])
    outfile.write("}")

