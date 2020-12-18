#%matplotlib inline
##!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import math as mt
import os

dirname = os.path.dirname(os.path.realpath(__file__))

def TimetoDetection(sigrate,sigrate_upper,sigrate_lower,backrate,maxdays,plot,save,name):

	"""
	Calculcate the time to detect CNO solar neutrinos to 3 and 5 \u03C3.
	Inputs needed: daily signal rate, daily background rate, maximum number of days to look for, plot = Yes or No, save plot = Yes or No, name to save plot as (don't include file type e.g. .pdf)
	Outputs: Time to 3 and 5 \u03C3, plot of significance of signal (if Yes selected)
	"""

	dirname = os.path.dirname(os.path.realpath(__file__))
	
	mus1 = sigrate
	mus1_upper = sigrate_upper
	mus1_lower = sigrate_lower
	mubg = backrate
	maxdays = maxdays
	sumsums1 = np.zeros(maxdays)
	sumsums1_upper = np.zeros(maxdays)
	sumsums1_lower = np.zeros(maxdays)
	sumsumbg = np.zeros(maxdays)
	loop = 1000
	
	for ntry in range(loop):
		if ntry % 10 == 0:
			print("Completion:", int(ntry/loop * 100),"%") #Print completion status of loop
		s1 = np.random.poisson(mus1,maxdays) #Poisson distribution of signal
		s1_upper = np.random.poisson(mus1_upper,maxdays)
		s1_lower = np.random.poisson(mus1_lower,maxdays)
		bg = np.random.poisson(mubg,maxdays) #Poisson distribution of background
		day = np.arange(0,maxdays,1)
		sums1 = np.cumsum(s1) #Cumulatively sum the signal
		sums1_upper = np.cumsum(s1_upper)
		sums1_lower = np.cumsum(s1_lower)
		sumbg = np.cumsum(bg) #Cumulatively sum the background
		sumsums1 = sumsums1 + sums1
		sumsums1_upper = sumsums1_upper + sums1_upper
		sumsums1_lower = sumsums1_lower + sums1_lower
		sumsumbg = sumsumbg + sumbg
    
	as1 = sumsums1/1000.0
	as1_upper = sumsums1_upper/1000.0
	as1_lower = sumsums1_lower/1000.0
	abg = sumsumbg/1000.0

	#plt.hist(abg,bins=maxdays,label='background')
	#plt.hist(as1,bins=maxdays,label='signal')
	#plt.legend()
	#plt.show()
	#plt.hist(abg,label='background')
	#plt.legend()
	#plt.show()

	print("Determining signal distribution")
	r = as1/(np.sqrt(abg))
	r_upper = as1_upper/(np.sqrt(abg))
	r_lower = as1_lower/(np.sqrt(abg))

	print("Finding position of 3\u03C3")
	horiz_line_data_3 = np.array([3 for i in range(len(day))]) #Horizontal line at 3 sigma
	print("Finding position of 5\u03C3")
	horiz_line_data_5 = np.array([5 for i in range(len(day))]) #Horizontal line at 5 sigma
	
	day_point_3 = np.arange(0,3.6,3/5)
	day_point = np.arange(0,6,1)

	print("Finding time to 3\u03C3")
	time_3 = np.argwhere(np.diff(np.sign(r - horiz_line_data_3))).flatten() #Time to 3 sig in days
	time_3_upper = np.argwhere(np.diff(np.sign(r_upper - horiz_line_data_3))).flatten()
	time_3_lower = np.argwhere(np.diff(np.sign(r_lower - horiz_line_data_3))).flatten()
	print("Finding time to 5\u03C3")
	time_5 = np.argwhere(np.diff(np.sign(r - horiz_line_data_5))).flatten() #Time to 5 sig in days
	time_5_upper = np.argwhere(np.diff(np.sign(r_upper - horiz_line_data_5))).flatten()
	time_5_lower = np.argwhere(np.diff(np.sign(r_lower - horiz_line_data_5))).flatten()

	vert_line_data_3 = np.array([time_3 for i in range(len(day_point))]) #Vertical line from where signal reach 3 sigma
	vert_line_data_3_upper = np.array([time_3_upper for i in range(len(day_point))])
	vert_line_data_3_lower = np.array([time_3_lower for i in range(len(day_point))])
	vert_line_data_5 = np.array([time_5 for i in range(len(day_point))]) #Vertical line from where signal reach 5 sigma
	vert_line_data_5_upper = np.array([time_5_upper for i in range(len(day_point))])
	vert_line_data_5_lower = np.array([time_5_lower for i in range(len(day_point))])
	
	if not time_3_lower and not time_5_lower: #Significance doesn't reach 3 sigma (or 5 sigma)

		print("The data takes > ", maxdays, "days to reach 3\u03C3")

		if plot == 'Yes' or plot == 'yes' or plot == 'Y' or plot == 'y':
			print("Plotting data")
			plt.plot(day,r,'r',label='CNO Signal (Mean Flux)')
			plt.plot(day,r_upper,'k',label='CNO Signal (GS98)')
			plt.plot(day,r_lower,'y',label='CNO Signal (AGSS09)')

			print("Formatting plot")
			plt.legend(loc='upper left')
			plt.xlabel('Days')
			plt.ylabel('Significance [\u03C3]')
			plt.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
			plt.ylim(0,1.3*max(r))
			plt.title('Time to see CNO Solar Neutrino Signal')
			
			if save == 'Yes' or save == 'yes' or save == 'Y' or save == 'y':
				print('Saving figure as', name,'.pdf')
				plt.savefig(dirname + '/../results/' + name + '.pdf')

			else:
				print('Not saving figure')
				
			plt.show()

		else:
			print('Plot not selected')

		if save == 'Yes' or save == 'yes' or save == 'Y' or save == 'y':		
			print('Saving data as', name,'.txt')
			with open(dirname + '/../results/' + name + '.txt', 'w+') as outfile:
				outfile.write("Signal Rate = %.3f per day \n" % sigrate)
				outfile.write("Signal Rate (High Metallicty) = %.3f per day \n" % sigrate_upper)
				outfile.write("Signal Rate (AGSS09) = %.3f per day \n" % sigrate_lower)
				outfile.write("Background Rate = %.3f per day \n" % backrate)
				outfile.write("Signal didn't reach 3\u03C3 in %i days" % maxdays)

		else:
			print('Not saving data')
		
	elif not time_5_lower: #Significance doesn't reach 5 sigma
		
		if plot == 'Yes' or plot == 'yes' or plot == 'Y' or plot == 'y':

			print("Plotting data")
			plt.plot(day,r,'r',label='CNO Signal')
			plt.plot(day,r_upper,'k',label='CNO Signal (GS98)')
			plt.plot(day,r_lower,'y',label='CNO Signal (AGSS09)')
			print("Plotting 3\u03C3")
			plt.plot(day, horiz_line_data_3, 'b--',label='3\u03C3 ("Evidence")')
			plt.plot(vert_line_data_3, day_point_3, 'b--') #'r--',label='3\u03C3 (Mean Flux)')
			plt.plot(vert_line_data_3_upper, day_point_3, 'b--') #'k--',label='3\u03C3 (GS98)')
			plt.plot(vert_line_data_3_lower, day_point_3,  'b--')#'y--',label='3\u03C3 (AGSS09)')

			print("Formatting plot")
			plt.legend(loc='upper left')
			plt.xlabel('Days')
			plt.ylabel('Significance [\u03C3]')
			plt.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
			plt.xlim(0,1.3*time_3_lower)
			plt.ylim(0,5)#1.3*max(r))
			plt.title('Time to see CNO Solar Neutrino Signal')

			if save == 'Yes' or save == 'yes' or save == 'Y' or save == 'y':
				
				print('Saving Figure as', name,'.pdf')
				plt.savefig(dirname + '/../results/' + name + '.pdf')

			else:
				
				print('Not saving figure')

			plt.show()

		else:

			print('Plot not selected')

		if save == 'Yes' or save == 'yes' or save == 'Y' or save == 'y':		
			print('Saving data as', name,'.txt')
			with open(dirname + '/../results/' + name + '.txt', 'w+') as outfile:
				outfile.write("Signal Rate = %.3f per day \n" % sigrate)
				outfile.write("Signal Rate (High Metallicty) = %.3f per day \n" % sigrate_upper)
				outfile.write("Signal Rate (AGSS09) = %.3f per day \n" % sigrate_lower)
				outfile.write("Background Rate = %.3f per day \n" % backrate)
				outfile.write("Time to 3\u03C3 (Mean Flux) = %i days \n" % time_3)
				outfile.write("Time to 3\u03C3 (Mean Flux) = %.3f years \n" % (time_3/365.25))
				outfile.write("Time to 3\u03C3 (GS98) = %i days \n" % time_3_upper)
				outfile.write("Time to 3\u03C3 (GS98) = %.3f years \n" % (time_3_upper/365.25))
				outfile.write("Time to 3\u03C3 (AGSS09) = %i days \n" % time_3_lower)
				outfile.write("Time to 3\u03C3 (AGSS09) = %.3f years \n" % (time_3_lower/365.25))
				outfile.write("Signal did't reach 5\u03C3 in %i days" % maxdays)

		else:
			print('Not saving data')

		print ("Time to 3\u03C3 (Mean Flux):", time_3, "days")
		print ("Time to 3\u03C3 (Mean Flux):", time_3/(365.25),"years")
		print ("Time to 3\u03C3 (GS98):", time_3_upper, "days")
		print ("Time to 3\u03C3 (GS98):", time_3_upper/(365.25),"years")
		print ("Time to 3\u03C3 (AGSS09):", time_3_lower, "days")
		print ("Time to 3\u03C3 (AGSS09):", time_3_lower/(365.25),"years")
		print("The data takes > ", maxdays, "days to reach 5\u03C3")
	

	elif max(r_lower)>vert_line_data_3_lower.all() and max(r_lower)>vert_line_data_5_lower.all(): #Significance passes 5 sigma
	
		if plot == 'Yes' or plot == 'yes' or plot == 'Y' or plot == 'y':

			print("Plotting data")
			plt.plot(day,r,'r',label='CNO Signal')
			plt.plot(day,r_upper,'k',label='CNO Signal (GS98)')
			plt.plot(day,r_lower,'y',label='CNO Signal (AGSS09)')
			print("Plotting 3\u03C3")
			plt.plot(day, horiz_line_data_3, 'b--',label='3\u03C3 ("Evidence")')
			plt.plot(vert_line_data_3, day_point_3,  'b--')#'r--',label='3\u03C3 (Mean Flux)')
			plt.plot(vert_line_data_3_upper, day_point_3,  'b--')#'k--',label='3\u03C3 (GS98)')
			plt.plot(vert_line_data_3_lower, day_point_3, 'b--') #'y--',label='3\u03C3 (AGSS09)')

			print("Plotting 5\u03C3")
			plt.plot(day, horiz_line_data_5, 'g--',label='5\u03C3 ("Discovery")')
			plt.plot(vert_line_data_5, day_point, 'g--')#'r--',label='5\u03C3 (Mean Flux)')
			plt.plot(vert_line_data_5_upper, day_point, 'g--')#'k--',label='5\u03C3 (GS98)')
			plt.plot(vert_line_data_5_lower, day_point, 'g--')#'y--',label='5\u03C3 (AGSS09)')

			print("Formatting plot")
			plt.legend(loc='upper left')
			plt.xlabel('Days')
			plt.ylabel('Significance [\u03C3]')
			plt.ticklabel_format(style='sci',axis='x', scilimits=(0,0))
			plt.xlim(0,1.1*time_5_lower)
			plt.ylim(0,10)#1.3*max(r))
			plt.title('Time to see CNO Solar Neutrino Signal')

			if save == 'Yes' or save == 'yes' or save == 'Y' or save == 'y':
				
				print('Saving Figure as', name,'.pdf')
				plt.savefig(dirname + '/../results/' + name + '.pdf')

			else:

				print('Not saving figure')

			plt.show()

		else:
			
			print('Plot not selected')

		if save == 'Yes' or save == 'yes' or save == 'Y' or save == 'y':		
			print('Saving data as', name,'.txt')
			with open(dirname + '/../results/' + name + '.txt', 'w+') as outfile:
				outfile.write("Signal Rate = %.3f per day \n" % sigrate)
				outfile.write("Signal Rate (High Metallicty) = %.3f per day \n" % sigrate_upper)
				outfile.write("Signal Rate (AGSS09) = %.3f per day \n" % sigrate_lower)
				outfile.write("Background Rate = %.3f per day \n" % backrate)
				outfile.write("Time to 3\u03C3 (Mean Flux) = %i days \n" % time_3)
				outfile.write("Time to 3\u03C3 (Mean Flux) = %.3f years \n" % (time_3/365.25))
				outfile.write("Time to 3\u03C3 (GS98) = %i days \n" % time_3_upper)
				outfile.write("Time to 3\u03C3 (GS98) = %.3f years \n" % (time_3_upper/365.25))
				outfile.write("Time to 3\u03C3 (AGSS09) = %i days \n" % time_3_lower)
				outfile.write("Time to 3\u03C3 (AGSS09) = %.3f years \n" % (time_3_lower/365.25))
				outfile.write("Time to 5\u03C3 (Mean Flux) = %i days \n" % time_5)
				outfile.write("Time to 5\u03C3 (Mean Flux) = %.3f years \n" % (time_5/365.25))
				outfile.write("Time to 5\u03C3 (GS98) = %i days \n" % time_5_upper)
				outfile.write("Time to 5\u03C3 (GS98) = %.3f years \n" % (time_5_upper/365.25))
				outfile.write("Time to 5\u03C3 (AGSS09) = %i days \n" % time_5_lower)
				outfile.write("Time to 5\u03C3 (AGSS09) = %.3f years \n" % (time_5_lower/365.25))

		else:
			print('Not saving data')

		print ("Time to 3\u03C3 (Mean Flux):", time_3, "days")
		print ("Time to 3\u03C3 (Mean Flux):", time_3/(365.25),"years")
		print ("Time to 3\u03C3 (GS98):", time_3_upper, "days")
		print ("Time to 3\u03C3 (GS98):", time_3_upper/(365.25),"years")
		print ("Time to 3\u03C3 (AGSS09):", time_3_lower, "days")
		print ("Time to 3\u03C3 (AGSS09):", time_3_lower/(365.25),"years")
		print ("Time to 5\u03C3 (Mean Flux):", time_5, "days")
		print ("Time to 5\u03C3 (Mean Flux):", time_5/(365.25),"years")
		print ("Time to 5\u03C3 (GS98):", time_5_upper, "days")
		print ("Time to 5\u03C3 (GS98):", time_5_upper/(365.25),"years")
		print ("Time to 5\u03C3 (AGSS09):", time_5_lower, "days")
		print ("Time to 5\u03C3 (AGSS09):", time_5_lower/(365.25),"years")

	else:
		print("Something has gone wrong determining time")

	return time_3,time_5

sigrate = float(eval(input('Daily signal rate (mean flux): ')))
sigrate_upper = float(eval(input('Daily signal rate (GS98): ')))
sigrate_lower = float(eval(input('Daily signal rate (AGSS09): ')))
backrate = float(eval(input('Daily background rate: ')))
maxdays = int(eval(input('Maximum number of days: ')))
plot = input('Plot signal? (Yes/No) ' )
save = input('Save figure? (Yes/No) ' )
name = input('Figure name (exclude file type, leave blank if not saving figure): ')

time_3,time_5=TimetoDetection(sigrate,sigrate_upper,sigrate_lower,backrate,maxdays,plot,save,name)
