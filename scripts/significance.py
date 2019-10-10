#%matplotlib inline

import matplotlib.pyplot as plt
import numpy as np
import math as mt

# create the data 
mus1 = 16.35 # signal rate
mubg = 250000 # background rate
sumsums1 = np.zeros(9000) # 100 experiment days - you should increase this for CNO solar neutrinos!
sumsumbg = np.zeros(9000)

for ntry in range(1000):
    s1 = np.random.poisson(mus1,9000)
    bg = np.random.poisson(mubg,9000)
    day = np.arange(0,9000,1)
    sums1 = np.cumsum(s1)
    sumbg = np.cumsum(bg)
    sumsums1 = sumsums1 + sums1
    sumsumbg = sumsumbg + sumbg
    
as1 = sumsums1/1000.0
abg = sumsumbg/1000.0

r = as1/(np.sqrt(as1+abg))

time = 8500

horiz_line_data = np.array([3 for i in range(len(day))])
day_point = np.arange(0,6,1)
vert_line_data = np.array([time for i in range(len(day_point))])
plt.plot(vert_line_data, day_point, 'r--')
plt.plot(day, horiz_line_data, 'r--') 
plt.plot(day,r,'r',label='CNO solar with Combined Watchman Backgrounds') # edit as required
plt.legend(loc='best')
plt.xlabel('days')
plt.ylabel('significance')
plt.title('Time to see CNO solar signal')
plt.show()


print ("Time:", time, "days")
print ("Time:", time/(365),"years")
