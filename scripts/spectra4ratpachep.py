import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_hep = 5.28e3


# read data for the hep chain neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
hep = np.genfromtxt(dirname + "/../data/hep.dat", delimiter="  ") # sns.ias.edu

Ehep = hep[:, 0]
fhep = hep[:, 1]*fi_hep

with open(dirname + "/../data/hep.ratdb",'w') as outfile:
    outfile.write("{\n")
    outfile.write('name: "hep"\n')
    outfile.write("valid_begin: [0,0],\n")
    outfile.write("valid_end:[0,0],\n")
    outfile.write("emin: %.5g\n" % Ehep[0])
    outfile.write("emax: %.5g\n" % Ehep[-1])
    outfile.write("spec_e: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%.5g," % Ehep[i])
        else: outfile.write("%.5g],\n" % Ehep[i])
    outfile.write("spec_flux: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%.5g," % fhep[i])
        else: outfile.write("%.5g],\n" % fhep[i])
    outfile.write("}")

