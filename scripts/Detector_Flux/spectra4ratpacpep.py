import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_pep = 7.98e8

Epep = 1.44 # Solar Neutrinos Spectroscopy 1704.06331

#EPEP=[Epep]*240
fPEP=np.linspace(0,fi_pep,240) #Create 240 points for hep spectrum between min and max flux at fixed energy

EPEP=np.linspace(Epep-0.5*(240/10000000),Epep+0.5*(240/10000000),240)

Rfid = (10026.35-3080)/10 #Radius fiducial volume in cm
hfid = 2*Rfid #Height of fiducial volume in cm

theta = 1.3103 * np.pi/180 #Average Solar elevation at Boubly

Aeff = (np.cos(theta) * hfid * 2 * Rfid) + (np.sin(theta) * np.pi * Rfid**2) #Effective area of fiducial volume exposed to Sun

#fPEP = np.pi * Rfid**2 * fPEP #Flux in s^-1 MeV^-1 in detector

fPEP = Aeff * fPEP #Flux in s^-1 MeV^-1 in detector

FreeProtons = 0.668559 # 0.668559 * 10^32 Free protons per kton of water
nktons = 2 * np.pi * (hfid/100) * (Rfid/100)**2 /1000 # Number of ktons of water in detector. 1000 m^3 in 1 kton
TNU = FreeProtons * nktons #Using s^-1 not year^-1 as this is what watchmakers uses.

FPEP = fPEP * TNU * EPEP #Flux in terms of TNU (using s^-1)

n_bins = len(EPEP)

#write the combined spectrum to a file in the required format for ratpac
with open(dirname + "/../../data/Detector_Flux/pep.ratdb",'w') as outfile:
    outfile.write("{\n")
    outfile.write('name: "SPECTRUM",\n')
    outfile.write('index: "pep",\n')
    outfile.write("valid_begin: [0,0],\n")
    outfile.write("valid_end:[0,0],\n")
    outfile.write("emin: %f\n" % Epep)
    outfile.write("emax: %f\n" % Epep)
    outfile.write("spec_e: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % EPEP[i])
        else: outfile.write("%f],\n" % EPEP[i])
    outfile.write("spec_mag: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % FPEP[i])
        else: outfile.write("%f],\n" % FPEP[i])
    outfile.write("}")

