#%matplotlib inline

import matplotlib.pyplot as plt
import numpy as np
import math as mt

# create the data 
mus1 = 16.35 # signal rate
mubg = 250000 # background rate

def TimetoDetection(sigrate,backrate):
	mus1 = sigrate
	musbg = backrate
	sumsums1 = np.zeros(100000)
	sumsumbg = np.zeros(100000)

	for ntry in range(1000):
		s1 = np.random.poisson(mus1,100000)
		bg = np.random.poisson(mubg,100000)
		day = np.arange(0,100000,1)
		sums1 = np.cumsum(s1)
		sumbg = np.cumsum(bg)
		sumsums1 = sumsums1 + sums1
		sumsumbg = sumsumbg + sumbg
    
	as1 = sumsums1/1000.0
	abg = sumsumbg/1000.0

	r = as1/(np.sqrt(as1+abg))

#time = 22000

	horiz_line_data_3 = np.array([3 for i in range(len(day))])
	horiz_line_data_5 = np.array([5 for i in range(len(day))])
	
	day_point = np.arange(0,6,1)
	
	time_3 = np.argwhere(np.diff(np.sign(r - horiz_line_data_3))).flatten() #Time to 3 sig in days
	time_5 = np.argwhere(np.diff(np.sign(r - horiz_line_data_5))).flatten() #Time to 5 sig in days
	vert_line_data_3 = np.array([time_3 for i in range(len(day_point))])
	vert_line_data_5 = np.array([time_5 for i in range(len(day_point))])
	
	plt.plot(vert_line_data_3, day_point, 'b--',label='3\u03C3')
	plt.plot(day, horiz_line_data_3, 'b--')
	plt.plot(vert_line_data_5, day_point, 'g--',label='5\u03C3')
	plt.plot(day, horiz_line_data_5, 'g--')  
	plt.plot(day,r,'r',label='CNO solar with Combined Watchman Backgrounds') # edit as required

	plt.legend(loc='upper left')
	plt.xlabel('days')
	plt.ylabel('significance')
	plt.title('Time to see CNO solar signal')
	plt.show()
	
	print ("Time to 3\u03C3:", time_3, "days")
	print ("Time to 3\u03C3:", time_3/(365),"years")
	print ("Time to 5\u03C3:", time_5, "days")
	print ("Time to 5\u03C3:", time_5/(365),"years")

	return

TimetoDetection(mus1,mubg)
