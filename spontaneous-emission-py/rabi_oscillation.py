from condition import *
import numpy as np

#  Rabi oscillation modeling.
k = 0  # k should be in (0.1, 0.5)
omegaWave = (1 - k) * omega21  # Hz (red light)
delta = omegaWave - omega21  # = k * omega21

omegaRabi = mu12 * elField0 / hBar
print(delta)
omega = 0.5 * (delta**2 + omegaRabi**2)**0.5

rhoTimeDepend11 = []
rhoTimeDepend22 = []
timeStep = timeStop / pntNumber

for mT in range(pntNumber):
    t = mT * timeStep
    rho22 = (0.5 * omegaRabi**2 * (1 - np.cos(2 * omega * t))) / (delta**2 + omegaRabi**2)
    rho11 = 1 - rho22
    print(rho11, ' ', rho22, ';')
    rhoTimeDepend22.append(rho22)
    rhoTimeDepend11.append(rho11)

print("Success!")
