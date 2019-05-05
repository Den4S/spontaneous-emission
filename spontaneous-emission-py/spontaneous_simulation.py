import random
import numpy as np
from matplotlib import pyplot as plt
eps = 0.01 #thresold value of photon radiation
N = 50000 # number of atoms
t = [] # array for radiation time T1, T2, ..., TN
gamma_spont = 1 # constant coefficient of spontaneous emission
global_time = 0
func_array = []
for p in np.arange(1,N+1):
    step = 2
    count = 0
    while step > eps: # finding moment of radiation
        step=random.random()
        count = count + 1
    M = (N -2*(p-1))/2
    gamma_current = gamma_spont*(N/2+M)*(N/2-M+1)
    global_time = global_time + (eps/gamma_current)*count
    func_array.append(gamma_current*global_time)
    t.append(global_time)

plt.plot(t,func_array)
plt.axvline(1/N*np.log(N),0,1)
plt.plot()
plt.show()

