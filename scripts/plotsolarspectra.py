import numpy as np
import matplotlib.pyplot as plt
import pdb

#integrated fluxes from N Vingoles
fi_pp = 5.98e10
fi_be1 = 5e8 #~~~~~placeholder value~~~~~
fi_be2 = 1.44e9
fi_pep = 7.98e8
fi_n = 4.93e8
fi_o = 5.46e8
fi_b = 2.78e6
fi_f = 2.05e6
fi_hep = 5.28e3

#read data for the nitrogen-13 neutrinos. Flux in cm^-2 s^-1 at Earth's surface, E in MeV
n13 = np.genfromtxt("n13.dat", delimiter="  ") # sns.ias.edu
o15 = np.genfromtxt("o15.dat", delimiter="  ") # ''
f17 = np.genfromtxt("f17.dat", delimiter="  ") # ''
hep = np.genfromtxt("hep.dat", delimiter="  ") # ''
b8 = np.genfromtxt("b8.dat", usecols=[0, 1]) # ''
pp_raw = np.genfromtxt("pp.dat") # ''
be7_raw = np.genfromtxt("be7.dat") # ''
En = n13[:, 0]
fn = n13[:, 1]*fi_n
Eo = o15[:, 0]
fo = o15[:, 1]*fi_o
Ef = f17[:, 0]
ff = f17[:, 1]*fi_f
Ehep = hep[:, 0]
fhep = hep[:, 1]
Eb = b8[:, 0]
fb = b8[:, 1]
Ebe2 = (be7_raw[:92, 0] / 1000) + .8613 #converting from keV to MeV then adding 861.3 keV to return to energy
fbe2 = be7_raw[:92, 1]*fi_be2
Ebe1 = (be7_raw[92:, 0] / 1000) + .3843
fbe1 = be7_raw[92:, 1]*fi_be1


#Ebe1 = 0.3843 # sns.ias.edu
#Ebe2 = 0.8613 # ''
Epep = 1.44 # Solar Neutrinos Spectroscopy 1704.06331

col1 = np.array(np.concatenate((pp_raw[:, 0], pp_raw[:, 2], pp_raw[:, 4], pp_raw[:, 6])))
col2 = np.array(np.concatenate((pp_raw[:, 1], pp_raw[:, 3], pp_raw[:, 5], pp_raw[:, 7])))
ppn = np.column_stack((col1, col2))


Epp = ppn[:, 0]
fpp = ppn[:, 1]

#En = np.concatenate((En,np.linspace(En[-1]+.05,Ef[-1],300)))
En = np.concatenate((En,Eo[200:]))
fn = np.concatenate((fn,np.linspace(0,0,300)))
#pdb.set_trace()

#group into energy bins
def IQR(dist):
    return np.percentile(dist,75) - np.percentile(dist,25)
#n_bins = int(np.ceil((max(En[-1], Eo[-1], Ef[-1])-min(En[0], Eo[0], Ef[0]))/(2*IQR(Eo)*len(Eo)**(-1/3)))) #Freedman-Diaconis
#n_bins = int(round(len(Eo)**(1/2)))
n_bins = 240
#pdb.set_trace()
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
    #pdb.set_trace()


Enbinned = Ebin - ((max(En[-1], Eo[-1], Ef[-1]) - min(En[0], Eo[0], Ef[0])) / n_bins)
Enbinned = Enbinned[1:]

fig,ax = plt.subplots()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1, 21)
ax.set_ylim(10, 10e12)
ax.plot(Enbinned, np.array(fbinned), 'kx', label='CNO')

#ax.plot(Eo, fo/Eo*fi_o, label='O', linestyle='-')
#ax.plot(En, fn/En'''*fi_n''', label='N')
#ax.plot(Ef, ff/Ef*'''fi_f''', label='F', linestyle=':')
#ax.plot(Ehep, fhep/Ehep'''*fi_hep''', label='hep')
#ax.plot(Eb, fb/Eb'''*fi_b''', label='8B')
#ax.plot(Epp, fpp/Epp'''*fi_pp''', label='pp')

ax.plot(Eo, fo, label='O', linestyle='--')
ax.plot(En, fn, label='N')
ax.plot(Ef, ff, label='F', linestyle=':')
ax.plot(Ehep, fhep*fi_hep, label='hep')
ax.plot(Eb, fb*fi_b, label='8B')
ax.plot(Epp, fpp*fi_pp, label='pp')
ax.plot(Ebe1, fbe1, label = '7Be')
ax.plot(Ebe2, fbe2, label = '7Be')
#ax.axvline(Ebe1, label='7Be (1)')
#ax.axvline(Ebe2, label='7Be (2)')
#ax.axvline(Epep, ymax = np.log(fi_pep), label='pep')
ax.plot([Epep,Epep],[0,fi_pep*0.5],label = 'pep')
ax.set_xlabel("Energy (MeV)")
ax.set_ylabel("Flux (cm^-2 s^-1)")
ax.legend()
ax.plot()
#ax.plot(Eo,fo,'b--')
#ax.plot(Ef,ff,'r-.')

plt.show()
