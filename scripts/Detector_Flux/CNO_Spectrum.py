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

import warnings
warnings.filterwarnings("ignore")

dirname = os.path.dirname(os.path.realpath(__file__))


#########################################################################################################


#integrated fluxes from N Vingoles
fi_n = 2.04e8 #2.78e8 #2.41e8#4.93e8
fi_o = 1.44e8 #2.05e8 #1.745e8#5.46e8
fi_f = 3.26e6 #5.92e8 #4.59e6#2.05e6


#########################################################################################################



def CNO_data(fi_n,fi_o,fi_f):

	"""
	Read "data" for the nitrogen-13, oxygen-15 and fluorine-17 neutrinos. Flux in cm^-2 s^-1 MeV^-1 at Earth's surface, E in MeV
	"Data" from http://www.sns.ias.edu/~jnb/SNdata/sndata.html

	Inputs: fi_n (Integrated N-13 flux in cm^-2 s^-1), fi_o (Integrated O-15 flux in cm^-2 s^-1), fi_f (Integrated F-17 flux in cm^-2 s^-1)

	Outputs: En (N-13 energy, MeV), fn (N-13 flux, cm^-2 s^-1 MeV^-1), Eo (O-15 energy, MeV), fo (O-15 flux, cm^-2 s^-1 MeV^-1), Ef (F-17 energy, MeV), ff (F-17 flux, cm^-2 s^-1 MeV^-1), n_bins

	"""


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

	return(En,fn,Eo,fo,Ef,ff,n_bins)



#########################################################################################################



def CNOfunc(CNO_data,fi_n,fi_o,fi_f):#(En,fn,Eo,fo,Ef,ff,n_bins):

	'''
	Combine N, O, F components to form CNO spectrum
	Inputs: Energy + flux of nitrogen, oxygen and fluorine neutrinos in order, number of bins for energy
	Output: Energy and flux of CNO neutrinos
	'''

	En,fn,Eo,fo,Ef,ff,n_bins = CNO_data(fi_n,fi_o,fi_f)

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


#########################################################################################################



def CNOPlot(CNOfunc,En,fn,Eo,fo,Ef,ff,n_bins):

	'''
	Plot CNO spectrum + Components
	Inputs: Function to combine CNO components, energy + flux of nitrogen, oxygen and fluorine neutrinos in order, number of bins for energy
	Output: Plot of CNO spectrum
	'''

	CNO=CNOfunc(En,fn,Eo,fo,Ef,ff,n_bins)
	#print("Max fn = %.4e " % max(fn))

	plt.plot(CNO[0], CNO[1], label='CNO')
	plt.plot(En,fn,label='$^{13}$N')
	plt.plot(Eo,fo,label='$^{15}$O')
	plt.plot(Ef,ff,label='$^{17}$F')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(0.1, 21)
	plt.ylim(10, 10e12)
	plt.xlabel("Energy [MeV]")
	#plt.ylabel("Flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]")
	plt.ylabel("Flux [cm$^{-2}$ s$^{-1}$]")
	plt.legend()
	plt.title('CNO Solar Neutrino Spectrum')
	#plt.savefig(dirname + "/../plots/CNOSolarNeutrinoSpectrum.pdf")
	plt.show()

	return

#CNOPlot(CNOfunc,En,fn,Eo,fo,Ef,ff,n_bins)



#########################################################################################################



Enbinned, fbinned=CNOfunc(CNO_data,fi_n,fi_o,fi_f)

FreeProtons = 0.668559 # 0.668559 * 10^32 Free protons per kton of water
FreeElectrons = 3.34e32
nktons=1.32 #Approximate fiducial volume
TNU = FreeProtons * nktons #Using s^-1 not year^-1 as this is what watchmakers uses.

#sigma_0 = 8.806e-45 #cm**2 https://scholarworks.umass.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1109&context=dissertations_2

Fbinned_int = integrate.trapz(fbinned,Enbinned)

fupper = 0.523e9 #CNO neutrino grand prix figure 8, table 1.3 of Direct measurements of the pp solar neutrinos in Borexino
flower = 0.337e9
fmean2 = (fupper+flower)/2

print("Calculated flux = %.5e cm^-2 s^-1" % Fbinned_int)
print("High metallicity flux = %.5e cm^-2 s^-1" % fupper)
print("Mean flux = %.5e cm^-2 s^-1" % fmean2)
print("Low metallicity flux = %.5e cm^-2 s^-1" % flower)



#########################################################################################################


def cross_section_ES(E):

	"""
	Calculate elastic scattering cross-section of electron neutrinos and electrons in cm**2
	Use of equation 1 in MINERvA, Measurement of neutrino flux from neutrino-electron elastic scattering, Physical Review D, Vol 93, 112007, 2016
	Method recommended by Teppei Katori, UCL

	Input: Energy (MeV) [array]
	Output: Cross-section of elastic scattering (cm**2) [array]

	dsigma/dy = Gf**2 * s/pi * [Cll**2 + Clr**2(1-y)**2]
	y = Te/Ev -> assume Te << Ev at room temperature
	Cll = Clr = 1/2 + sin**2(theta_w) for nu_e
	Cll = 1/2 - sin**2(theta_w) for nu_mu, nu_tau
	Clr = sin**2(theta_w) for nu_mu, nu_tau
	Assume s = Ev * me

	"""


	me = 0.511e-3 #GeV
	Ev = E/1000 #GeV
	Gf = 1.1663787e-5 #GeV
	sin2theta_w = 0.23

	dsigma = Gf**2 * me * Ev * 2 * (0.5 + sin2theta_w)**2 / np.pi
	sigma_cm = 0.389e-27 * dsigma #convert to cm^2 from natural units
	sigma_err = np.std(sigma_cm)

	plt.plot(E,sigma_cm)
	plt.xlabel('Energy (MeV)')
	plt.ylabel('\u03C3 (cm$^{2}$)')
	plt.title('Cross-Section of Elastic Scattering of \u03BD$_e$e$^-$')
	#plt.savefig(dirname + '/../../results/Cross-Section_enu.pdf')
	plt.show()

	return sigma_cm,sigma_err


########################################################################################################


def energy_hist_load(CNO_data,fi_n,fi_o,fi_f):

	energy_hist = np.genfromtxt(dirname + "/particle_energy_values_CNO.txt", delimiter=",")

	print(len(energy_hist))

	En,fn,Eo,fo,Ef,ff,n_bins = CNO_data(fi_n,fi_o,fi_f)

	(n,bins,patches)=plt.hist(energy_hist,n_bins)
	plt.xlabel('Energy (MeV)')
	plt.ylabel('Events')
	plt.show()

	return(energy_hist,n,bins,patches)



#########################################################################################################
energy_hist,n,bins,patches=energy_hist_load(CNO_data,fi_n,fi_o,fi_f)
En,fn,Eo,fo,Ef,ff,n_bins = CNO_data(fi_n,fi_o,fi_f)
sigma_cm,sigma_err = cross_section_ES(Enbinned)
print("SIG_ERR =", sigma_err)


def CNO_rate(sigma_cm,E,energy_hist,n_bins,fupper,fmean,flower):


	year_rate = 3600*24*365.25
	FreeProtons = 0.668559 # 0.668559 * 10^32 Free protons per kton of water
	FreeElectrons = 3.34e32
	nktons=1.32 #Approximate fiducial volume


	rate_upper = 2*year_rate*FreeElectrons*sigma_cm*fupper*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
	rate_mean = 2*year_rate*FreeElectrons*sigma_cm*fmean2*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
	rate_lower = 2*year_rate*FreeElectrons*sigma_cm*flower*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins

	fitted_upper=gaussian_filter1d(rate_upper, sigma=2)
	rate_int_upper = integrate.trapz(rate_upper,Enbinned)

	fitted_mean=gaussian_filter1d(rate_mean, sigma=2)
	rate_int_mean = integrate.trapz(rate_mean,Enbinned)

	fitted_lower=gaussian_filter1d(rate_lower, sigma=2)
	rate_int_lower = integrate.trapz(rate_lower,Enbinned)

	rate_model = 1.5*year_rate*FreeElectrons*sigma_cm*fbinned*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
	rate_model_int = integrate.trapz(rate_model,Enbinned)

	print("Rate (high metallicity) = %.4e /s" % (rate_int_upper/0.02/365.25/3600/24))
	print("Rate (high metallicity) = %.4e /day" % (rate_int_upper/0.02/365.25))
	print("Rate (mean flux) = %.4e /s" % (rate_int_mean/0.02/365.25/3600/24))
	print("Rate (mean flux) = %.4e /day" % (rate_int_mean/0.02/365.25))
	print("Rate (low metallicity) = %.4e /s" % (rate_int_lower/0.02/365.25/3600/24))
	print("Rate (low metallicity) = %.4e /day" % (rate_int_lower/0.02/365.25))

	rate_upper_25 = rate_upper*25
	rate_mean_25 = rate_mean*25
	rate_lower_25 = rate_lower*25

	fitted_upper_25=gaussian_filter1d(rate_upper_25, sigma=2)
	fitted_mean_25=gaussian_filter1d(rate_mean_25, sigma=2)
	fitted_lower_25=gaussian_filter1d(rate_lower_25, sigma=2)

	return


def CNO_rate_plot():

	fig,ax = plt.subplots()
	ax.plot(Enbinned,fitted_upper,label='CNO (GS98)')
	ax.plot(Enbinned,fitted_mean,label='CNO')
	ax.plot(Enbinned,fitted_lower,label='CNO (AGSS09)')
	plt.xlim(0.6,2*max(Enbinned))
	plt.ylim(min(fitted_upper),5e4)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Energy (MeV)')
	plt.ylabel('Events / 0.02 MeV / kton / year')
	ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
	ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
	ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
	ax.xaxis.set_minor_formatter(StrMethodFormatter('{x:.1f}'))
	plt.legend()
	plt.title('Expected Energy Spectra')
	#plt.savefig(dirname + '/../../results/Signal_per_kton_CNO.pdf')
	plt.show()


	fig,ax = plt.subplots()
	ax.plot(Enbinned,fitted_upper_25,label='CNO (GS98)')
	ax.plot(Enbinned,fitted_mean_25,label='CNO')
	ax.plot(Enbinned,fitted_lower_25,label='CNO (AGSS09)')
	plt.xlim(0.6,2*max(Enbinned))
	plt.ylim(min(fitted_upper_25),5e4)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Energy (MeV)')
	plt.ylabel('Events / 0.02 MeV / year')
	ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
	ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
	ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
	ax.xaxis.set_minor_formatter(StrMethodFormatter('{x:.1f}'))
	plt.legend()
	plt.title('Expected Energy Spectra for 25 kiloton Detector')
	#plt.savefig(dirname + '/../../results/Signal_25kton_CNO.pdf')
	plt.show()


	return


#rate_plot_upper = 2*3600*24*365.25*FreeElectrons*sigma_cm*fupper*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
#rate_plot_mean = 2*3600*24*365.25*FreeElectrons*sigma_cm*fmean2*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
#rate_plot_lower = 2*3600*24*365.25*FreeElectrons*sigma_cm*flower*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
#rate_plot_unit=rate_plot*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins

#fitted_upper=gaussian_filter1d(rate_plot_upper, sigma=2)
#rate_plot_int_upper = integrate.trapz(rate_plot_upper,Enbinned)

#fitted_mean=gaussian_filter1d(rate_plot_mean, sigma=2)
#rate_plot_int_mean = integrate.trapz(rate_plot_mean,Enbinned)

#fitted_lower=gaussian_filter1d(rate_plot_lower, sigma=2)
#rate_plot_int_lower = integrate.trapz(rate_plot_lower,Enbinned)

#print('Int rate = ', (rate_plot_int)/0.02/365.25/3600/24)

#rate_plot_model = 1.5*3600*24*365.25*FreeElectrons*sigma_cm*fbinned*n[0]/(0.02*Enbinned)/len(energy_hist)/n_bins
#rate_plot_model_int = integrate.trapz(rate_plot_model,Enbinned)

#print("Rate (high metallicity) = %.4e /s" % (rate_plot_int_upper/0.02/365.25/3600/24))
#print("Rate (high metallicity) = %.4e /day" % (rate_plot_int_upper/0.02/365.25))
#print("Rate (mean flux) = %.4e /s" % (rate_plot_int_mean/0.02/365.25/3600/24))
#print("Rate (mean flux) = %.4e /day" % (rate_plot_int_mean/0.02/365.25))
#print("Rate (low metallicity) = %.4e /s" % (rate_plot_int_lower/0.02/365.25/3600/24))
#print("Rate (low metallicity) = %.4e /day" % (rate_plot_int_lower/0.02/365.25))

#fig,ax = plt.subplots()
#ax.plot(Enbinned,fitted_upper,label='CNO (GS98)')
#ax.plot(Enbinned,fitted_mean,label='CNO')
#ax.plot(Enbinned,fitted_lower,label='CNO (AGSS09)')
#plt.xlim(0.6,2*max(Enbinned))
#plt.ylim(min(fitted_upper),5e4)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('Energy (MeV)')
#plt.ylabel('Events / 0.02 MeV / kton / year')
#ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
#ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
#ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
#ax.xaxis.set_minor_formatter(StrMethodFormatter('{x:.1f}'))
#plt.legend()
#plt.title('Expected Energy Spectra')
#plt.savefig(dirname + '/../../results/Signal_per_kton_CNO.pdf')
#plt.show()




#########################################################################################################




#rate_plot_upper_25 = rate_plot_upper*25
#rate_plot_mean_25 = rate_plot_mean*25
#rate_plot_lower_25 = rate_plot_lower*25

#fitted_upper_25=gaussian_filter1d(rate_plot_upper_25, sigma=2)
#fitted_mean_25=gaussian_filter1d(rate_plot_mean_25, sigma=2)
#fitted_lower_25=gaussian_filter1d(rate_plot_lower_25, sigma=2)

#fitted=gaussian_filter1d(rate_plot_unit, sigma=2)

#fig,ax = plt.subplots()
#ax.plot(Enbinned,fitted_upper_25,label='CNO (GS98)')
#ax.plot(Enbinned,fitted_mean_25,label='CNO')
#ax.plot(Enbinned,fitted_lower_25,label='CNO (AGSS09)')
#plt.xlim(0.6,2*max(Enbinned))
#plt.ylim(min(fitted_upper_25),5e4)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('Energy (MeV)')
#plt.ylabel('Events / 0.02 MeV / year')
#ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
#ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
#ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
#ax.xaxis.set_minor_formatter(StrMethodFormatter('{x:.1f}'))
#plt.legend()
#plt.title('Expected Energy Spectra for 25 kiloton Detector')
#plt.savefig(dirname + '/../../results/Signal_25kton_CNO.pdf')
#plt.show()



#########################################################################################################



#plt.plot(Enbinned,sigma_cm)
#plt.xlabel('Energy (MeV)')
#plt.ylabel('\u03C3 (cm$^{2}$)')
#plt.title('Cross-Section of Elastic Scattering of \u03BD$_e$e$^-$')
#plt.savefig(dirname + '/../../results/Cross-Section_enu.pdf')
#plt.show()




#########################################################################################################




x = Enbinned
y = fbinned/400#/400#*500#/n_bins

ysmoothed = gaussian_filter1d(y, sigma=2)

y_int = integrate.trapz(ysmoothed,Enbinned)
print("Calculated flux = %.5e cm^-2 s^-1" % y_int)

#ysmoothed[0].append(0)
ysmoothed = np.append(ysmoothed, 0)
x = np.append(x,max(Ef))

fig,ax = plt.subplots()
#plt.plot(Enbinned,fbinned,label='CNO')
ax.plot(En,fn,label='$^{13}$N')
ax.plot(Eo,fo,label='$^{15}$O')
ax.plot(Ef,ff,label='$^{17}$F')
ax.plot(x, ysmoothed,label='CNO')
plt.xlabel('Energy (MeV)')
plt.ylabel('Flux (cm$^{-2}$ s$^{-1}$ MeV$^{-1}$)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.6,2*max(Enbinned))
#plt.ylim(1e6,1e10)
ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
ax.xaxis.set_minor_formatter(StrMethodFormatter('{x:.1f}'))
#plt.grid()
plt.legend()
#plt.savefig(dirname + '/../../results/Expected_Flux.pdf')
plt.show()



########################################################################################################



fbinned=2*fbinned*(fmean2/integrate.trapz(fbinned/500,Enbinned))/500#*4.2e8/(integrate.trapz(fbinned,Enbinned)*400)
print("int fbinned %.4e" % (integrate.trapz(fbinned,Enbinned)))

rate_model_mean = FreeElectrons*sigma_cm*(fbinned+0.66e8)/2#ysmoothed[0:-1]/2 The factor of 1/2 is for oscillations
rate_model_int_mean = integrate.trapz(rate_model_mean,Enbinned)
rate_model_upper = FreeElectrons*sigma_cm*((fupper+0.73e8)/fmean2)*fbinned/2#ysmoothed[0:-1]/2
rate_model_int_upper = integrate.trapz(rate_model_upper,Enbinned)
rate_model_lower = FreeElectrons*sigma_cm*(flower+0.52e8)/fmean2*fbinned/2#ysmoothed[0:-1]/2
rate_model_int_lower = integrate.trapz(rate_model_lower,Enbinned)

print("Rate (high metallicity) = %.4e /s" % (rate_model_int_upper))
print("Rate (high metallicity) = %.4e /day" % (rate_model_int_upper*3600*24))
print("Rate (mean flux) = %.4e /s" % (rate_model_int_mean))
print("Rate (mean flux) = %.4e /day" % (rate_model_int_mean*3600*24))
print("Rate (low metallicity) = %.4e /s" % (rate_model_int_lower))
print("Rate (low metallicity) = %.4e /day" % (rate_model_int_lower*3600*24))

#u_mean = rate_model_int_mean * np.sqrt(np.var(rate_model_mean) + sum(FreeElectrons*sigma_cm*fbinned/2)/rate_model_int_mean**2)



#######################################################################################################




Fel_err = 0.03e32
fbinned_err = 0.66e8#np.std(fbinned)
fbinned_err_high = 0.73e8
fbinned_err_low = 0.53e8

#print(Fel_err**2/FreeElectrons**2)
#print(fbinned_err**2/fbinned**2)
#print(sigma_err**2/sigma_cm**2)

u = rate_model_int_mean * np.sqrt(np.var(rate_model_mean) * (Fel_err * fbinned_err * sigma_err)**2 / rate_model_int_mean**2) #http://www.m4ssl.npl.co.uk/wp-content/uploads/2012/02/Determining-the-uncertainty-associated-with-integrals-of-spectral-quantities.pdf equation 7.6
print(u*3600*24)






###################################################################################################






#with open("SolarNeutrinoRates.txt",'a+') as outfile:
#	outfile.write("Rate for CNO in 1.32 kTon fiducial volume = %e per second\n" % ratemean)
#	outfile.write("Rate for CNO (low metallicity) in 1.32 kTon fiducial volume = %e per second\n" % ratelower)
#	outfile.write("Rate for CNO (high metallicity) in 1.32 kTon fiducial volume = %e per second\n" % rateupper)

#write the combined spectrum to a file in the required format for ratpac
#with open(dirname + "/../../data/Detector_Flux/CNO.ratdb",'w') as outfile:
#with open(dirname + "/Test_CNO_Data.ratdb",'w') as outfile:
#    outfile.write("{\n")
#    outfile.write('name: "SPECTRUM",\n')
#    outfile.write('index: "CNO",\n')
#    outfile.write("valid_begin: [0,0],\n")
#    outfile.write("valid_end:[0,0],\n")
#    outfile.write("emin: %f\n" % Enbinned[0])
#    outfile.write("emax: %f\n" % Enbinned[-1])
#    outfile.write("spec_e: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % Enbinned[i])
#        else: outfile.write("%f],\n" % Enbinned[i])
#    outfile.write("spec_mag: [")
#    for i in range(n_bins-1):
#        if i != (n_bins-2): outfile.write("%f," % fbinned[i])
#        else: outfile.write("%f],\n" % fbinned[i])
#    outfile.write("}")

