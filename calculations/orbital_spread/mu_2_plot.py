# Imports:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import pandas as pd
import scienceplots


# Style for plotting:
plt.style.use(['science','ieee'])


# Read data
data = pd.read_csv("mu_vals.txt", sep=";")


# Remove whitespace
data.columns = data.columns.str.strip()


# Get results
# cc-pVDZ
distances = data["R"]                                                        # Å 


# mus from cl-CCSD calculations
mu_2_ccpvdz     = data["mu_2_ccpvdz"]                                       
mu_2_ccpvtz = data["mu_2_ccpvtz"]                                    
mu_2_aug_o_only = data["mu_2_aug_o_only"]                                    


# Plotting
fig, ax = plt.subplots()

plt.xlabel("Distance [Å]")
plt.ylabel(r"2nd order orbital spread [Å]")

ax.plot(distances, mu_2_ccpvdz, label = "cc-pVDZ cl-CCSD", color = "purple", linestyle='-', linewidth=0.2, ms = 2, marker = '^', mfc='none', mew=0.3, markevery=2)                     # cc-pVDZ             CCSD
ax.plot(distances, mu_2_ccpvdz, label = "cc-pVTZ cl-CCSD", color = "green", linestyle='-', linewidth=0.2, ms = 2, marker = 'o', mfc='none', mew=0.3, markevery=2)                     # cc-pVDZ             CCSD
ax.plot(distances, mu_2_aug_o_only, label = "aug-cc-pVTZ / cc-pVTZ cl-CCSD", color = "orange", linestyle='-', linewidth=0.2, ms = 2, marker = 's',  mfc='none', mew=0.3, markevery=2)              # cc-pVDZ             CL-CCSD


ax.ticklabel_format(axis='y', style='plain', useOffset=False)

ax.legend(ncol=2, fontsize=4.5, frameon=False, loc='upper right')

plt.savefig("orbital_spread.pdf", dpi = 600)
