pntNumber = 10000  # number of atoms

# time variable for modeling
timeStop = 10  # s


eCharge = 1.60217662 * 10**(-19)  # charge of electron
angsr = 10**(-10)  # angstrom in meters
hBar = 1.054571800 * 10**(-34)  # Plank constant

elField0 = 1  # amplitude value of external field

energy1 = -13.6 * eCharge
energy2 = energy1 / 4  # Bohr model energy
omega21 = (energy1 - energy2) / hBar

mu12 = angsr * eCharge  # magnetic moment
mu21 = mu12

