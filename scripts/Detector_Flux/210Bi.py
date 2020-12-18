import numpy as np
import matplotlib.pyplot as plt
import pdb
import os
from scipy import integrate
from scipy import optimize
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from scipy.ndimage.filters import gaussian_filter1d
import seaborn as sns

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
fbinned=fbinned/500

energy_hist_210Bi = np.genfromtxt(dirname + "/particle_energy_values_210Bi.txt", delimiter=",")

#print((energy_hist_210Bi))

energy_hist_CNO = np.genfromtxt(dirname + "/particle_energy_values_CNO_bonsai.txt", delimiter=",")

#print(len(energy_hist_CNO))


energy_hist_210Bi = energy_hist_210Bi# * 2500 / len(energy_hist_210Bi)
(n_210Bi,bins_210Bi,patches_210Bi)=plt.hist(energy_hist_210Bi,n_bins)
plt.xlabel('Energy (MeV)')
plt.ylabel('Events')
plt.show()

#print(n_210Bi)

E = np.linspace(0,1.17,n_bins)

#print(n[0])

(n_CNO,bins_CNO,patches_CNO)=plt.hist(energy_hist_CNO,n_bins)
plt.xlabel('Energy (MeV)')
plt.ylabel('Events')
#plt.xlim(0.6,1.4)
plt.show()

Bi = n_210Bi[0]# *2500/len(energy_hist_210Bi)
CNO = n_CNO[0]# * 1e-3/len(energy_hist_CNO)

plt.plot(E,n_CNO[0])
plt.plot(E,n_210Bi[0])
plt.show()

plt.plot(E,Bi)
plt.plot(E,CNO)
plt.show()

energy_hist = np.genfromtxt(dirname + "/particle_energy_values_CNO.txt", delimiter=",")

(n,bins,patches)=plt.hist(energy_hist,n_bins)
plt.xlabel('Energy (MeV)')
plt.ylabel('Events')
#plt.xlim(0.6,1.4)
plt.show()

FreeElectrons = 3.32e32
fmean = 4.2e8
fupper = 0.523e9 #CNO neutrino grand prix figure 8, table 1.3 of Direct measurements of the pp solar neutrinos in Borexino
flower = 0.337e9
fmean2 = (fupper+flower)/2
me = 0.511e-3 #GeV
Ev = Enbinned/1000 #GeV
Gf = 1.1663787e-5 #GeV
sin2theta_w = 0.23

dsigma = Gf**2 * me * Ev * 2 * (0.5 + sin2theta_w)**2 / np.pi #Paper sent by Teppei
sigma_cm = 0.389e-27 * dsigma #convert to cm^2 from natural units

rate_plot_upper = 2*3600*24*365.25*FreeElectrons*sigma_cm*fupper*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
rate_plot_mean = 2*3600*24*365.25*FreeElectrons*sigma_cm*fmean2*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
rate_plot_lower = 2*3600*24*365.25*FreeElectrons*sigma_cm*flower*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
#rate_plot_unit=rate_plot*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins

fitted_upper=gaussian_filter1d(rate_plot_upper, sigma=2)
rate_plot_int_upper = integrate.trapz(rate_plot_upper,Enbinned)

fitted_mean=gaussian_filter1d(rate_plot_mean, sigma=2)
rate_plot_int_mean = integrate.trapz(rate_plot_mean,Enbinned)

fitted_lower=gaussian_filter1d(rate_plot_lower, sigma=2)
rate_plot_int_lower = integrate.trapz(rate_plot_lower,Enbinned)

Bi = Bi#*2500*3600*24*365.25/0.02/E/len(energy_hist_210Bi)/n_bins

fig,ax = plt.subplots()
ax.plot(Enbinned,fitted_upper,label='CNO (High Metallicity)')
ax.plot(Enbinned,fitted_mean,label='CNO')
ax.plot(Enbinned,fitted_lower,label='CNO (Low Metallicity)')
plt.plot(E,Bi,label='210Bi')
plt.xlim(0.6,2*max(Enbinned))
plt.ylim(min(fitted_upper),5e4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (MeV)')
plt.ylabel('Events / 0.02 Mev / kton / year')
ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
ax.xaxis.set_minor_formatter(StrMethodFormatter('{x:.1f}'))
plt.legend()
plt.title('Expected energy spectra in water Cherenkov detector')
#plt.savefig(dirname + '/../../results/Signal_per_kton_CNO.pdf')
plt.show()
