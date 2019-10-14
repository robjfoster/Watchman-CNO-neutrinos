import numpy as np
import matplotlib.pyplot as plt
import pdb
import os

dirname = os.path.dirname(os.path.realpath(__file__))

#integrated fluxes from N Vingoles
fi_n = 4.93e8
fi_o = 5.46e8
fi_f = 2.05e6


# read data for the nitrogen-13 neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
n13 = np.genfromtxt(dirname + "/../data/n13.dat", delimiter="  ") # sns.ias.edu
o15 = np.genfromtxt(dirname + "/../data/o15.dat", delimiter="  ") # ''
f17 = np.genfromtxt(dirname + "/../data/f17.dat", delimiter="  ") # ''
comspec = np.genfromtxt(dirname+"/../data/CNO_old.ratdb", delimiter=",")
En = n13[:, 0]
fn = n13[:, 1]#*fi_n
Eo = o15[:, 0]
fo = o15[:, 1]#*fi_o
Ef = f17[:, 0]
ff = f17[:, 1]#*fi_f

Ec = comspec[0,:]
fc = comspec[1,:]


En = np.concatenate((En, Eo[200:]))
fn = np.concatenate((fn, np.linspace(0, 0, 300)))

# group into energy bins
def IQR(dist):
    return np.percentile(dist,75) - np.percentile(dist,25)
# n_bins = int(np.ceil((max(En[-1], Eo[-1], Ef[-1])-min(En[0], Eo[0], Ef[0]))/(2*IQR(Eo)*len(Eo)**(-1/3)))) #Freedman-Diaconis
# n_bins = int(round(len(Eo)**(1/2))) # sqrt N
n_bins = 240
Ebin = np.linspace(min(En[0], Eo[0], Ef[0]), max(En[-1], Eo[-1], Ef[-1]), n_bins)
fbinned = []
En_bin = np.linspace(En[0], En[-1], n_bins)
fn_binned = []
for i in range(n_bins - 1):
    fbin = 0
    for E, f in zip(En, fn):
        if Ebin[i] < E < Ebin[i+1]:
            fbin += f
    for E, f in zip(Eo, fo):
        if Ebin[i] < E < Ebin[i+1]:
            fbin += f
    for E, f in zip(Ef, ff):
        if Ebin[i] < E < Ebin[i+1]:
            fbin += f
    fbinned.append(fbin)

Enbinned = Ebin - ((max(En[-1], Eo[-1], Ef[-1]) - min(En[0], Eo[0], Ef[0])) / n_bins)
Enbinned = Enbinned[1:]

fig,ax = plt.subplots()
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlim(0.1, 21)
#ax.set_ylim(10, 10e12)
ax.plot(Enbinned, np.array(fbinned), 'kx', label='New')
ax.plot(Ec,fc,'ro',label='Summer')

ax.set_xlabel("Energy (MeV)")
ax.set_ylabel("Flux (cm^-2 s^-1)")
ax.legend()
ax.plot()
#plt.show()
plt.savefig(dirname + "/../plots/comparesolarspectra.png")
