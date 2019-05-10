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


eps = 0.01  # thresold value of photon radiation
N = 10000  # number of atoms
radT = []  # array for radiation time T1, T2, ..., TN
gammaSpont = 1  # constant coefficient of spontaneous emission
globalTime = 0
funcArray = []
for p in np.arange(1, N+1):
    step = 2
    count = 0
    while step > eps:  # finding moment of radiation
        step = random.random()
        count = count + 1
    M = (N - 2 * (p - 1)) / 2
    gammaCurrent = gammaSpont * (N / 2 + M) * (N / 2 - M + 1)
    globalTime = globalTime + (eps/gammaCurrent)*count
    funcArray.append(gammaCurrent * globalTime)
    radT.append(globalTime)

funcMax = max(funcArray)  # funcArray extremum
indMax = funcArray.index(funcMax)

# let's find: 2 * deltaT is full width at half maximum
flagL = False
flagR = False
indL = 0
indR = 0
for i in range(N + 1):
    if (not flagL) and (funcArray[i] > funcMax / 2):
        indL = i
        flagL = True
    if (not flagR) and flagL and (funcArray[i] < funcMax / 2):
        indR = i
        flagR = True
deltaT = (radT[indR] - radT[indL]) / 2
print("deltaT =", deltaT, "s")  # a parameter for Lorentz function

fig, ax = plt.subplots()  # creating a plot
ax.set_xlabel('t, s', size=12)
ax.set_ylabel('$\gamma T_n$', size=12)

plt.plot(radT, funcArray, 'blue')  # our function
plt.plot(radT, lorentz(radT, radT[indMax], deltaT, max(funcArray)), 'green', linestyle='-', linewidth=1)  # Lorentz

plt.plot([1 / N * np.log(N), 1 / N * np.log(N)], [0, 1.05 * funcMax], 'magenta',
         linestyle='--', linewidth=1)  # theoretical extremum line
plt.plot([radT[indMax], radT[indMax]], [0, 1.05 * funcMax], 'red', linestyle='-.', linewidth=1)  # extremum line

plt.legend(("$\gamma T_n(t)$", "Lorentz", "1 / (N * ln(N))", "$\max \; (\gamma T_n(t))$"), loc=1)  # plot legend

plt.plot()
plt.show()

