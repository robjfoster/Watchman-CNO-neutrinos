import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_be2 = 1.44e9

# read data for the second be-7 set of neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
be2 = np.genfromtxt(dirname + "/../data/be7_2.dat", delimiter="  ") # sns.ias.edu
Ebe2 = (be2[:, 0]/1000) + .8613 # KeV to MeV and addition to return to energy
fbe2 = be2[:, 1]*fi_be2

n_bins = len(Ebe2)

#write the combined spectrum to a file in the required format for ratpac
with open(dirname + "/../data/be2.ratdb",'w') as outfile:
    outfile.write("{\n")
    outfile.write('name: "SPECTRUM",\n')
    outfile.write('index: "7Be_2",\n')
    outfile.write("valid_begin: [0,0],\n")
    outfile.write("valid_end:[0,0],\n")
    outfile.write("emin: %f\n" % Ebe2[0])
    outfile.write("emax: %f\n" % Ebe2[-1])
    outfile.write("spec_e: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % Ebe2[i])
        else: outfile.write("%f],\n" % Ebe2[i])
    outfile.write("spec_mag: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % fbe2[i])
        else: outfile.write("%f],\n" % fbe2[i])
    outfile.write("}")

