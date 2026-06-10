# Import 
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# Read data
data = pd.read_csv("coupling_elements_01_0-1_data_aug_DZ.txt", sep=";")


# Remove whitespace
data.columns = data.columns.str.strip()


# Get results
# cc-pVDZ
distances   = data["R"]                 # Å 

H_AB      = data["H_AB"]                # eV  


# exp fitting 
def exponential(x, a, b, c):
    return a * np.exp(-b * x) + c


popt, pcov = curve_fit(exponential, distances, H_AB, p0=(1, 0.1, 1))


# Plotting 
fig, ax = plt.subplots()

plt.xlabel("q [Å]")
plt.ylabel(r"$V$ [eV]")



ax.scatter(distances, H_AB, label=r"$|V_{01}^{0,-1}|$, aug-cc-pVDZ cl-CCSD-eom", color="red", marker='^', s=6, facecolors='none', linewidths=0.3)
ax.plot(distances, exponential(distances, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt), color = 'black', lw = 0.6)

plt.legend()
plt.show()
