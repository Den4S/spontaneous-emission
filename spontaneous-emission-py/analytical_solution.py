import random
import math
import numpy as np
from matplotlib import pyplot as plt


def findTau(time, function):  # let's find 2 * deltaT: it's full width at half maximum
    flagL = False
    flagR = False
    indL = 0
    indR = 0
    functionMax = max(function)
    for i in range(len(function)):
        if (not flagL) and (function[i] > functionMax / 2):
            indL = i
            flagL = True
        if (not flagR) and flagL and (function[i] < functionMax / 2):
            indR = i
            flagR = True
    delta = (time[indR] - time[indL]) / 2
    return delta


def MonteCarloMethod1(N):
    # eps - thresold value of photon radiation
    # N - number of atoms
    eps = 1
    radT = []  # array for radiation time T1, T2, ..., TN
    gammaSpont = 1  # constant coefficient of spontaneous emission
    globalTime = 0
    funcArray = []
    for p in np.arange(1, N + 1):
        step = 2
        count = 0
        while step > eps:  # finding moment of radiation
            step = random.random()
            count = count + 1
        M = (N - 2 * (p - 1)) / 2
        gammaCurrent = gammaSpont * (N / 2 + M) * (N / 2 - M + 1)
        globalTime = globalTime + (eps / gammaCurrent) * count
        funcArray.append(gammaCurrent * globalTime)
        radT.append(globalTime)

    funcMax = max(funcArray)  # funcArray extremum
    indMax = funcArray.index(funcMax)

    deltaT = findTau(radT, funcArray)

    return radT, funcArray, indMax, deltaT


def analytical(N):  # if $\epsilon == 1$ => analytical solution exists!
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
    tNphoton, tNintensity = analytical(N)  # T_N and I_N arrays for analytical solution
    tN, fN, iM, deltaMcT = MonteCarloMethod1(N)  # our Monte Carlo method for comparison with analytical solution

    normCoeff = fN[iM] / max(tNintensity)  # BAD normalization!!!

    # iMintens = tNintensity.index(max(tNintensity))  # artificial translation!!!
    # dT = tN[iM] - tNphoton[iMintens]

    for ind in range(len(tNintensity)):
        tNintensity[ind] = tNintensity[ind] * normCoeff
        # tNphoton[ind] = tNphoton[ind] + dT

    delta = findTau(tNphoton, tNintensity)

    return tNphoton, tNintensity, tN, fN, iM, delta, deltaMcT


# MAIN
nAtoms = 100000
tIntensArr, intensArr, tArr, funcArr, iMax, tau, tauMc = forEps1(nAtoms)

tauAnArr = []  # full width at half maximum for analytical solution
tauMcArr = []  # full width at half maximum for MC method
nAtArr = []
for nAtoms in range(1000, 100000, 1000):
    tIntensArrN, intensArrN, tArrN, funcArrN, iMaxN, tauN, tauMcN = forEps1(nAtoms)
    tauAnArr.append(tauN)
    tauMcArr.append(tauMcN)
    nAtArr.append(nAtoms)

# displaying a plot for one value of nAtoms
figAn, ax = plt.subplots(figsize=(8, 5))  # creating a plot
ax.set_xlabel('t, s', size=12)
ax.set_ylabel('$\gamma_{m}$', size=12)

plt.plot(tIntensArr, intensArr, 'm', linestyle='-', linewidth=1.5)  # intensity function (analytical)
plt.plot(tArr, funcArr, 'blue', linestyle='--', linewidth=1)  # our function
plt.plot([1 / nAtoms * np.log(nAtoms), 1 / nAtoms * np.log(nAtoms)], [0, 1.05 * max(funcArr)], 'm',
         linestyle='--', linewidth=1)  # theoretical extremum line
plt.plot([tArr[iMax], tArr[iMax]], [0, 1.05 * max(funcArr)], 'r', linestyle='-.', linewidth=1)  # extremum line

plt.legend(("$I_{analytical}$", "$I_{MC}$",
            "1 / (N * ln(N))", "$\max \; (\gamma T_n(t))$"), loc=1)  # plot legend

plt.title('$\gamma_{m}(T_{m})$')

plt.plot()
plt.show(block=False)

# plot tau(N) for analytical solution when
fig2, ax2 = plt.subplots(figsize=(8, 5))  # creating a plot
ax2.set_xlabel('Number of atoms', size=12)
ax2.set_ylabel('$\\tau$', size=12)

plt.plot(nAtArr, tauMcArr, 'b', linestyle='-', linewidth=1.5)  # tau(N) for MC method
plt.plot(nAtArr, tauAnArr, 'm', linestyle='--', linewidth=1.5)  # tau(N) for analytical solution
plt.title('$\\tau(N)$')
plt.legend(("$MC\; for\; \epsilon = 1$", "$Analytical\; for\; \epsilon = 1$"), loc=1)  # plot legend

plt.title('$\\tau(N)$')

plt.plot()
plt.show()
