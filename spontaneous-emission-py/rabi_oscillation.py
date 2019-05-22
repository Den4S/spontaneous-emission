from condition import *
import numpy as np
from matplotlib import pyplot as plt


#  Rabi oscillation modeling.
# k = 0.1  # k should be in (0.1, 0.5)
# omegaWave = (1 - k) * omega21  # Hz (red light)
# delta = omegaWave - omega21  # = k * omega21

omegaRabi = 1  # consider Rabi frequency = 1
delta = 0.2  # measure delta frecuency in Rabi frequences
omega = 0.5 * (delta**2 + omegaRabi**2)**0.5

rhoTimeDepend11 = []
rhoTimeDepend22 = []
population_inversion = []
time = []
timeStep = timeStop / pntNumber

for mT in range(pntNumber):
    t = mT * timeStep
    time.append(t)
    rho22 = (0.5 * omegaRabi**2 * (1 - np.cos(2 * omega * t))) / (delta**2 + omegaRabi**2)
    rho11 = 1 - rho22
    # print(rho11, ' ', rho22, ';')
    rhoTimeDepend22.append(rho22)
    rhoTimeDepend11.append(rho11)
    population_inversion.append(abs(rho11)**2-abs(rho22)**2)

fig = plt.figure()
ax1 = fig.add_subplot(121)
plt.plot(time, rhoTimeDepend22, 'blue')
plt.plot(time, rhoTimeDepend11, 'red')
plt.legend(("$c_{a_{1}}(t)$", "$c_{a_{2}}(t)$"), loc=1)
ax1.set_xlabel('time, a.u.', size=12)
ax1.set_ylabel('Probability amplitude', size=12)  # plotting probability amplitudes
ax1.set_title('Rabi oscillations ($\Delta$='+str(delta)+', $\Omega_{0}$)')
plt.grid()
ax2 = fig.add_subplot(122)
plt.plot(time, population_inversion, 'green')
ax2.set_xlabel('time, a.u.', size=12)
ax2.set_ylabel('$|c_{a{1}}|^{2}-|c_{a{2}}|^{2}$', size=12)
ax2.set_title('Population inversion')  # plotting population inversion/ time
# plt.legend('Population inversion')
plt.grid()
plt.show()

print("Success!")
