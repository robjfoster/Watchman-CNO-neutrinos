import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_pep = 7.98e8

Epep = 1.44 # Solar Neutrinos Spectroscopy 1704.06331

EPEP=[Epep]*240
fPEP=np.linspace(0,fi_pep,240) #Create 240 points for hep spectrum between min and max flux at fixed energy

n_bins = len(EPEP)

#write the combined spectrum to a file in the required format for ratpac
with open(dirname + "/../data/pep.ratdb",'w') as outfile:
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
        if i != (n_bins-2): outfile.write("%f," % fPEP[i])
        else: outfile.write("%f],\n" % fPEP[i])
    outfile.write("}")

