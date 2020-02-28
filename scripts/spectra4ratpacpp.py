import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles

fi_pp = 5.98e10


# read data for the pp chain neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
pp = np.genfromtxt(dirname + "/../data/pp.dat", delimiter="  ") # sns.ias.edu

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

n_bins = len(Epp)

#write the combined spectrum to a file in the required format for ratpac
with open(dirname + "/../data/pp.ratdb",'w') as outfile:
    outfile.write("{\n")
    outfile.write('name: "pp"\n')
    outfile.write("valid_begin: [0,0],\n")
    outfile.write("valid_Eppd:[0,0],\n")
    outfile.write("emin: %.5g\n" % Epp[0])
    outfile.write("emax: %.5g\n" % Epp[-1])
    outfile.write("spec_e: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%.5g," % Epp[i])
        else: outfile.write("%.5g],\n" % Epp[i])
    outfile.write("spec_flux: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%.5g," % fpp[i])
        else: outfile.write("%.5g],\n" % fpp[i])
    outfile.write("}")

