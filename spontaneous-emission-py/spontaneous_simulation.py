import random
import numpy as np
from matplotlib import pyplot as plt


def lorentz(x_arr, x_0, gamma, max_y):  # procedure for building a custom Lorentz function
    y_arr = []
    for x in x_arr:
        y_arr.append(gamma / (np.pi * ((x - x_0)**2 + gamma**2)))
    lorentz_max = max(y_arr)
    for m in range(len(y_arr)):
        y_arr[m] = y_arr[m] * max_y / lorentz_max  # Y scale
    return y_arr


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


def MonteCarloMethod(eps, N):
    # eps - thresold value of photon radiation
    # N - number of atoms
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

    # let's find: 2 * deltaT is full width at half maximum
    deltaT = findTau(radT, funcArray)

    lArr = lorentz(radT, radT[indMax], deltaT, max(funcArray))
    return radT, funcArray, lArr, indMax, deltaT


# MAIN
epsilon = 0.1
nAtomsFix = 10000

tArr, funcArr, lorentzArr, iMax, tau = MonteCarloMethod(epsilon, nAtomsFix)

tauArr = []
nAtArr = []
for nAtoms in range(1000, 100000, 1000):
    tArrN, funcArrN, lorentzArrN, iMaxN, tauN = MonteCarloMethod(epsilon, nAtoms)
    tauArr.append(tauN)
    nAtArr.append(nAtoms)

# displaying a plot $\gamma T_n(t)$ for one value of nAtoms
fig1, ax1 = plt.subplots()  # creating a plot
ax1.set_xlabel('t, s', size=12)
ax1.set_ylabel('$\gamma T_n$', size=12)

plt.plot(tArr, funcArr, 'blue')  # our function
plt.plot(tArr, lorentzArr, 'green', linestyle='-', linewidth=1)  # Lorentz

plt.plot([1 / nAtomsFix * np.log(nAtomsFix), 1 / nAtomsFix * np.log(nAtomsFix)], [0, 1.05 * max(funcArr)], 'magenta',
         linestyle='--', linewidth=1)  # theoretical extremum line
plt.plot([tArr[iMax], tArr[iMax]], [0, 1.05 * max(funcArr)], 'red', linestyle='-.', linewidth=1)  # extremum line

plt.legend(("$\gamma T_n(t)$", "Lorentz", "1 / (N * ln(N))", "$\max \; (\gamma T_n(t))$"), loc=1)  # plot legend

plt.plot()
plt.show(block=False)

# plot tau(N)
fig2, ax2 = plt.subplots()  # creating a plot
ax2.set_xlabel('Number of atoms', size=12)
ax2.set_ylabel('$tau$', size=12)

plt.plot(nAtArr, tauArr, 'red', linestyle='-.', linewidth=1.5)  # tau(N)

plt.legend(("$\epsilon =$" + str(epsilon),), loc=1)  # plot legend

plt.plot()
plt.show()
