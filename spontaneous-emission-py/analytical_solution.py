from condition import *
import math
from matplotlib import pyplot as plt


# if $\epsilon == 1$ => analytical solution exists!

tNphoton = []  # n photons emitted
tNintensity = []
gamma = 1

t0 = math.log(pntNumber) / (gamma * (pntNumber + 1))

for photonNumber in range(pntNumber-1):
    frstAtan = math.atanh((2 * (pntNumber / 2 - (photonNumber + 1)) - 1) / (pntNumber + 1))  # M = N / 2 - n
    scndAtan = math.atanh((2 * (pntNumber / 2) - 1) / (pntNumber + 1))  # M = N / 2
    tNew = -2 / (gamma * (pntNumber + 1)) * (frstAtan - scndAtan)
    tNphoton.append(tNew)
    intensityNew = gamma * ((pntNumber + 1) / 2)**2 * (1.0 / (math.cosh((pntNumber + 1) * (tNew - t0) / 2)))**2
    tNintensity.append(intensityNew)

fig, ax = plt.subplots()  # creating a plot
ax.set_xlabel('t, s', size=12)
ax.set_ylabel('$Intensity$', size=12)

plt.plot(tNphoton, tNintensity, 'm', linestyle='-', linewidth=1)  # intensity function
plt.legend(("$I_{analyt}(t)$",), loc=1)  # plot legend

plt.plot()
plt.show()

print(tNphoton)