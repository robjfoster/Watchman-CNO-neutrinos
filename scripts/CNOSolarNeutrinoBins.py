#Testing the effect of the number of bins used to aggregate CNO flux on the flux itself. Previous issue where more bins lowered flux solved by normalising flux

import numpy as np
import matplotlib.pyplot as plt
import pdb
import os
from scipy import integrate

dirname = os.path.dirname(os.path.realpath(__file__)) #File path

#integrated fluxes from N Vingoles
fi_n = 4.93e8
fi_o = 5.46e8
fi_f = 2.05e6

# read data for the nitrogen-13 neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
n13 = np.genfromtxt(dirname + "/../data/n13.dat", delimiter="  ") # sns.ias.edu
o15 = np.genfromtxt(dirname + "/../data/o15.dat", delimiter="  ") # ''
f17 = np.genfromtxt(dirname + "/../data/f17.dat", delimiter="  ") # ''

En = n13[:, 0]
fn = n13[:, 1] * fi_n #Multiply by integrated flux
Eo = o15[:, 0]
fo = o15[:, 1] * fi_o
Ef = f17[:, 0]
ff = f17[:, 1] * fi_f

En = np.concatenate((En, Eo[200:]))
fn = np.concatenate((fn, np.linspace(0, 0, 300)))

# group into energy bins
def IQR(dist):
    return np.percentile(dist,75) - np.percentile(dist,25)
n_bins = int(np.ceil((max(En[-1], Eo[-1], Ef[-1])-min(En[0], Eo[0], Ef[0]))/(2*IQR(Eo)*len(Eo)**(-1/3)))) 

def CNOfunc(En,fn,Eo,fo,Ef,ff,n_bins):

    Ebin = np.linspace(min(En[0], Eo[0], Ef[0]), max(En[-1], Eo[-1], Ef[-1]), n_bins)
    fbinned = []
    En_bin = np.linspace(En[0], En[-1], n_bins)
    fn_binned = []
    for i in range(n_bins - 1):
        fbin = 0
        for E, f in zip(En, fn):
            if Ebin[i] < E < Ebin[i+1]:
                fbin += f*n_bins # Need to multiply flux by number of bins to normalise
        for E, f in zip(Eo, fo):
            if Ebin[i] < E < Ebin[i+1]:
                fbin += f*n_bins
        for E, f in zip(Ef, ff):
            if Ebin[i] < E < Ebin[i+1]:
                fbin += f*n_bins
        fbinned.append(fbin)

    Enbinned = Ebin - ((max(En[-1], Eo[-1], Ef[-1]) - min(En[0], Eo[0], Ef[0])) / n_bins)
    Enbinned = Enbinned[1:]

    return Enbinned, fbinned

CNO=CNOfunc(En,fn,Eo,fo,Ef,ff,n_bins)


def CNOPlot(CNOfunc,En,fn,Eo,fo,Ef,ff,n_bins):

    CNO=CNOfunc(En,fn,Eo,fo,Ef,ff,n_bins)

    plt.plot(CNO[0], CNO[1], label=i)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1, 21)
    plt.ylim(10, 10e12)
    plt.xlabel("Energy [MeV]")
    plt.ylabel("Flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]")

n_bins_list=(5,n_bins,10,50,100,150,200,240,300)

for i in n_bins_list:
    CNOPlot(CNOfunc,En,fn,Eo,fo,Ef,ff,i) #Plot for several different bin sizes
plt.plot(CNO[0], CNO[1], 'kx', label="IQR Method (8)") #Compare to IQR method
plt.title('Number of Bins vs CNO Solar Neutrino Spectrum')
plt.legend(ncol=2)
plt.savefig(dirname + "/../plots/CNO_vs_n_bins_renomalised.pdf")
plt.show()
