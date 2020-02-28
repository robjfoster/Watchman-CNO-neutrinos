import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_b = 2.78e6

# read data for the boron-8 neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
b8 = np.genfromtxt(dirname + "/../data/b8.dat", delimiter="       ") # sns.ias.edu
Eb = b8[:, 0]
fb = b8[:, 1]*fi_b # b8.dat not delimited properly (replaced number of spaces in file to match)

#write the combined spectrum to a file in the required format for ratpac
with open(dirname + "/../data/b8.ratdb",'w') as outfile:
    outfile.write("{\n")
    outfile.write('name: "b8"\n')
    outfile.write("valid_begin: [0,0],\n")
    outfile.write("valid_end:[0,0],\n")
    outfile.write("emin: %.5g\n" % Eb[0])
    outfile.write("emax: %.5g\n" % Eb[-1])
    outfile.write("spec_e: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%.5g," % Eb[i])
        else: outfile.write("%.5g],\n" % Eb[i])
    outfile.write("spec_flux: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%.5g," % fb[i])
        else: outfile.write("%.5g],\n" % fb[i])
    outfile.write("}")

