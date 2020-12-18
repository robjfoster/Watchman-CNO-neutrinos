#import numpy as np
#import matplotlib.pyplot as plt
#import os

#dirname = os.path.dirname(os.path.realpath(__file__))

#energy_hist = np.genfromtxt(dirname + "/particle_energy_values.txt", delimiter=",")
#print(len(energy_hist))

#(n,bins,patches)=plt.hist(energy_hist,240)
#plt.xlim(0.6,1.4)

#print(np.shape(n))#,bins,len(bins))
#print((n[0]))

#plt.show()

import numpy as np
import matplotlib.pyplot as plt
import pdb
import os
from scipy import integrate
from scipy import optimize
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles
fi_n = 4.93e8
fi_o = 5.46e8
fi_f = 2.05e6

# read data for the nitrogen-13, oxygen-15 and fluorine-17 neutrinos. Flux in cm^-2 s^-1 MeV^-1 at Earth's surface, E in MeV
n13 = np.genfromtxt(dirname + "/../../data/n13.dat", delimiter="  ") # sns.ias.edu
o15 = np.genfromtxt(dirname + "/../../data/o15.dat", delimiter="  ") # ''
f17 = np.genfromtxt(dirname + "/../../data/f17.dat", delimiter="  ") # ''

En = n13[:, 0]
fn = n13[:, 1] * fi_n
Eo = o15[:, 0]
fo = o15[:, 1] * fi_o
Ef = f17[:, 0]
ff = f17[:, 1] * fi_f

n_bins=int(min(len(En),len(Eo),len(Ef))/2)

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
Enbinned, fbinned=CNOfunc(En,fn,Eo,fo,Ef,ff,n_bins)
fbinned = fbinned/500

plt.plot(Enbinned, fbinned)
plt.xscale('log')
plt.yscale('log')
plt.show()

me = 0.511e-3 #GeV
Ev = Enbinned/1000 #GeV
Gf = 1.1663787e-5 #GeV
sin2theta_w = 0.23

dsigma = Gf**2 * me * Ev * 2 * (0.5 + sin2theta_w)**2 / np.pi
sigma_cm = 0.389e-27 * dsigma
#sigma = integrate.trapz(dsigma,Ev)

plt.plot(Enbinned,sigma_cm)
plt.show()
