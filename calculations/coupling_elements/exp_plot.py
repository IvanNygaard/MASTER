# Imports:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import pandas as pd
import scienceplots


# Style for plotting:
plt.style.use(['science','ieee'])


# For reference
def exponential(H0, k, x):
    return H0 * np.exp(-k * x)


# Read data
data00_dz     = pd.read_csv("coupling_elements_00_0-1_data_DZ.txt", sep=";")
data01_dz     = pd.read_csv("coupling_elements_01_0-1_data_DZ.txt", sep=";")

data00_aug_dz = pd.read_csv("coupling_elements_00_0-1_data_aug_DZ.txt", sep=";")
data01_aug_dz = pd.read_csv("coupling_elements_01_0-1_data_aug_DZ.txt", sep=";")


# Remove whitespace
data00_dz.columns = data00_dz.columns.str.strip()
data01_dz.columns = data01_dz.columns.str.strip()

data00_aug_dz.columns = data00_aug_dz.columns.str.strip()
data01_aug_dz.columns = data01_aug_dz.columns.str.strip()


# Get results
# cc-pVDZ
distances   = data00_dz["R"]                                                        # Å 


# cc-pVDZ
H_AB_00      = data00_dz["H_AB"]                                             # eV  
H_AB_01_CISD = data01_dz["H_AB_CISD"] * 27.2114                              # eV (NB! Data given in au for CISD but eV for CCSD in the .txt file!)
H_AB_01_CCSD = data01_dz["H_AB"]                                             # eV


# Plotting
fig, ax = plt.subplots()

plt.xlabel("q [Å]")
plt.ylabel(r"$V$ [eV]")



ax.plot(distances, H_AB_01_CCSD, label = r"$|V_{01}^{0,-1}|$, cc-pVDZ cl-CCSD-eom", color = "red", linestyle='-', marker='^', linewidth=0.2, ms = 2, mfc='none', mew=0.3, markevery=2)               # aug-cc-pVDZ   cl-CCSD-eom
ax.plot(distances, H_AB_00, label = r"$|V_{00}^{0,-1}|$, cc-pVDZ cl-CCSD-eom", color = "blue", linestyle='-', marker='v', linewidth=0.2, ms = 2, mfc='none', mew=0.3, markevery=2)                   # aug-cc-pVDZ   cl-CCSD-eom
#ax.plot(distances, H_AB_01_CISD, label = r"$|V_{01}^{0,-1}|$, aug-cc-pVDZ CISD", color = "black", linestyle='-', marker='d', linewidth=0.2, ms = 2, mfc='none', mew=0.3, markevery = 2)                  # aug-cc-pVDZ   cl-CISD-eom
#ax.plot(distances, exponential(1.0, 1.65, distances), label = r"H$_0$exp(-$\beta$x)", color = "forestgreen", linestyle='--', marker='o', linewidth=0.2, ms = 2, mfc='none', mew=0.3, markevery=2)    # reference exp(-beta * x)

ax.ticklabel_format(axis='y', style='plain', useOffset=False)

ax.legend(ncol=2, fontsize=4.5, frameon=False)

plt.savefig("ccsd_v_cisd.pdf", dpi = 600)


# aug-cc-pVDZ
H_AB_00         = data00_aug_dz["H_AB"]     # eV  
H_AB_01         = data01_aug_dz["H_AB"]     # eV

fig, ax = plt.subplots()

plt.xlabel("q [Å]")
plt.ylabel(r"$V$ [eV]")


ax.plot(distances, H_AB_01, label = r"$|V_{01}^{0,-1}|$, aug-cc-pVDZ cl-CCSD-eom", color = "red", linestyle='-', marker='^', linewidth=0.2, ms = 2, mfc='none', mew=0.3, markevery=2)               # aug-cc-pVDZ   cl-CCSD-eom
ax.plot(distances, H_AB_00, label = r"$|V_{00}^{0,-1}|$, aug-cc-pVDZ cl-CCSD-eom", color = "blue", linestyle='-', marker='v', linewidth=0.2, ms = 2, mfc='none', mew=0.3, markevery=2)                   # aug-cc-pVDZ   cl-CCSD-eom
ax.plot(distances, H_AB_01_CISD, label = r"$|V_{01}^{0,-1}|$, aug-cc-pVDZ CISD", color = "black", linestyle='-', marker='d', linewidth=0.2, ms = 2, mfc='none', mew=0.3, markevery = 2)                  # aug-cc-pVDZ   cl-CISD-eom
ax.plot(distances, exponential(1.0, 1.65, distances), label = r"H$_0$exp(-$\beta$x)", color = "forestgreen", linestyle='--', marker='o', linewidth=0.2, ms = 2, mfc='none', mew=0.3, markevery=2)    # reference exp(-beta * x)

ax.ticklabel_format(axis='y', style='plain', useOffset=False)

ax.legend(ncol=2, fontsize=4.5, frameon=False)

plt.savefig("ccsd_aug_DZ_00_01.pdf", dpi = 600)
