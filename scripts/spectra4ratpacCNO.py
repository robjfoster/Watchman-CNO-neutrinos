import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

import numpy as np
import matplotlib.pyplot as plt
import pdb
import os
from scipy import integrate

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles
fi_n = 4.93e8
fi_o = 5.46e8
fi_f = 2.05e6


# read data for the nitrogen-13, oxygen-15 and fluorine-17 neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
n13 = np.genfromtxt(dirname + "/../data/n13.dat", delimiter="  ") # sns.ias.edu
o15 = np.genfromtxt(dirname + "/../data/o15.dat", delimiter="  ") # ''
f17 = np.genfromtxt(dirname + "/../data/f17.dat", delimiter="  ") # ''
comspec = np.genfromtxt(dirname+"/../data/CNO_old.ratdb", delimiter=",") #Summer students CNO spectrum
En = n13[:, 0]
fn = n13[:, 1] * fi_n
Eo = o15[:, 0]
fo = o15[:, 1] * fi_o
Ef = f17[:, 0]
ff = f17[:, 1] * fi_f

n_bins=int(min(len(En),len(Eo),len(Ef))/2)
#En = np.concatenate((En, Eo[200:]))
#fn = np.concatenate((fn, np.linspace(0, 0, 300)))

#n_bins=int(min(len(En),len(Eo),len(Ef))/2)

def CNOfunc(En,fn,Eo,fo,Ef,ff,n_bins):

    '''
    Combine N, O, F components to form CNO spectrum
    Inputs: Energy + flux of nitrogen, oxygen and fluorine neutrinos in order, number of bins for energy
    Output: Energy and flux of CNO neutrinos
    '''

    Ebin = np.linspace(min(En[0], Eo[0], Ef[0]), max(En[-1], Eo[-1], Ef[-1]), n_bins)
    fbinned = []
    En_bin = np.linspace(En[0], En[-1], n_bins)
    fn_binned = []
    for i in range(n_bins - 1):
        fbin = 0
        for E, f in zip(En, fn):
            if Ebin[i] < E < Ebin[i+1]:
                fbin += f*n_bins
        for E, f in zip(Eo, fo):
            if Ebin[i] < E < Ebin[i+1]:
                fbin += f*n_bins
        for E, f in zip(Ef, ff):
            if Ebin[i] < E < Ebin[i+1]:
                fbin += f*n_bins
        fbinned.append(fbin)
    Ebin = np.append(Ebin,max(Ebin))
    Enbinned = Ebin - ((max(En[-1], Eo[-1], Ef[-1]) - min(En[0], Eo[0], Ef[0])) / n_bins)
    Enbinned = Enbinned[1:]

    fbinned=np.append(fbinned,0)

    return Enbinned, fbinned

def CNOPlot(CNOfunc,En,fn,Eo,fo,Ef,ff,n_bins):

    '''
    Plot CNO spectrum + Components
    Inputs: Function to combine CNO components, energy + flux of nitrogen, oxygen and fluorine neutrinos in order, number of bins for energy
    Output: Plot of CNO spectrum
    '''

    CNO=CNOfunc(En,fn,Eo,fo,Ef,ff,n_bins)

    plt.plot(CNO[0], CNO[1], label='CNO')
    plt.plot(En,fn,label='$^{13}$N')
    plt.plot(Eo,fo,label='$^{15}$O')
    plt.plot(Ef,ff,label='$^{17}$F')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1, 21)
    plt.ylim(10, 10e12)
    plt.xlabel("Energy [MeV]")
    plt.ylabel("Flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]")
    plt.legend()
    plt.title('CNO Solar Neutrino Spectrum')
    plt.savefig(dirname + "/../plots/CNOSolarNeutrinoSpectrum.pdf")
    plt.show()

#CNOPlot(CNOfunc,En,fn,Eo,fo,Ef,ff,n_bins)

Enbinned, fbinned=CNOfunc(En,fn,Eo,fo,Ef,ff,n_bins)

#write the combined spectrum to a file in the required format for ratpac
with open(dirname + "/../data/CNO.ratdb",'w') as outfile:
    outfile.write("{\n")
    outfile.write('name: "SPECTRUM",\n')
    outfile.write('index: "CNO",\n')
    outfile.write("valid_begin: [0,0],\n")
    outfile.write("valid_end:[0,0],\n")
    outfile.write("emin: %f\n" % Enbinned[0])
    outfile.write("emax: %f\n" % Enbinned[-1])
    outfile.write("spec_e: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % Enbinned[i])
        else: outfile.write("%f],\n" % Enbinned[i])
    outfile.write("spec_mag: [")
    for i in range(n_bins-1):
        if i != (n_bins-2): outfile.write("%f," % fbinned[i])
        else: outfile.write("%f],\n" % fbinned[i])
    outfile.write("}")

