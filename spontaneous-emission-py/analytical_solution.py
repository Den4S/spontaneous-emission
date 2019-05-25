from condition import *
import math
import numpy as np
from matplotlib import pyplot as plt
import spontaneous_simulation as spsim


# if $\epsilon == 1$ => analytical solution exists!
def analytical(N):
    tNphoton = []  # n photons emitted
    tNintensity = []
    gamma = 1

    t0 = math.log(N) / (gamma * (N + 1))

    for photonNumber in range(N-1):
        frstAtanh = math.atanh((2 * (N / 2 - (photonNumber + 1)) - 1) / (N + 1))  # M = N / 2 - n
        scndAtanh = math.atanh((2 * (N / 2) - 1) / (N + 1))  # M = N / 2
        tNew = -2 / (gamma * (N + 1)) * (frstAtanh - scndAtanh)
        tNphoton.append(tNew)
        intensityNew = gamma * ((N + 1) / 2)**2 * (1.0 / (math.cosh((N + 1) * (tNew - t0) / 2)))**2
        tNintensity.append(intensityNew)

    return tNphoton, tNintensity


def forEps1(N):  # calculates and display plots for epsilon = 1
    eps = 1
    tNphoton, tNintensity = analytical(N)  # T_N and I_N arrays for analytical solution
    tN, fN, lorN, iM = spsim.MonteCarloMethod(eps, N)

    maxI = max(tNintensity) # BAD normalization!!!
    maxF = fN[iM]
    for ind in range(len(tNintensity)):
        tNintensity[ind] = tNintensity[ind] * maxF / maxI

    return tNphoton, tNintensity, tN, fN, lorN, iM


# MAIN
nAtoms = 100000
tN, iN, tArr, funcArr, lorentzArr, iMax = forEps1(nAtoms)

# display the plot
figAn, ax = plt.subplots()  # creating a plot
ax.set_xlabel('t, s', size=12)
ax.set_ylabel('$\gamma T_n$', size=12)

plt.plot(tN, iN, 'm', linestyle='-', linewidth=1)  # intensity function (analytical)
plt.plot(tArr, funcArr, 'blue')  # our function
plt.plot([1 / nAtoms * np.log(nAtoms), 1 / nAtoms * np.log(nAtoms)], [0, 1.05 * max(funcArr)], 'm',
         linestyle='--', linewidth=1)  # theoretical extremum line
plt.plot([tArr[iMax], tArr[iMax]], [0, 1.05 * max(funcArr)], 'red', linestyle='-.', linewidth=1)  # extremum line

plt.legend(("$I_{analitical}(t)$", "$\gamma T_n(t)$", "1 / (N * ln(N))",
            "$\max \; (\gamma T_n(t))$"), loc=1)  # plot legend

plt.plot()
plt.show()
