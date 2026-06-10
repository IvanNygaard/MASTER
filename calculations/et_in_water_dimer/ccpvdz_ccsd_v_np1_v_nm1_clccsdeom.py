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
cc_pVDZ_data = pd.read_csv("cc-pVDZ_CCSD_v_CL_CCSD_EOM_Np1_Nm1.txt", sep=";")


# Remove whitespace
cc_pVDZ_data.columns = cc_pVDZ_data.columns.str.strip()


# Get results
# cc-pVDZ
distances   = cc_pVDZ_data["R"]                                                        # Å 


# Energies (NOTE: Not interaction, just AB)
E_CCSD_AB      = cc_pVDZ_data["E_CCSD_AB"]  
E_CL_CCSD_AB   = cc_pVDZ_data["E_CL_CCSD"]
E_CCSD_Np1_AB  = cc_pVDZ_data["E_CL_CCSD_EOM_AB_Np1"]                                        
E_CCSD_Nm1_AB  = cc_pVDZ_data["E_CL_CCSD_EOM_AB_Nm1"]                                       


# Plotting
fig, ax = plt.subplots()

plt.xlabel("Distance [Å]")
plt.ylabel(r"Energy [$E_{\mathrm{h}}$]")


ax.plot(distances, E_CCSD_AB, label = "cc-pVDZ CCSD", color = "black", linestyle = '-', linewidth=0.2, ms = 2, marker = 'x',  mfc='none', mew=0.3, markevery=2)                                                          # cc-pVDZ             CCSD
ax.plot(distances, E_CL_CCSD_AB, label = "cc-pVDZ cl-CCSD", color = "red", linestyle = '-', linewidth=0.2, ms = 2, marker = 's',  mfc='none', mew=0.3, markevery=2) 

ax.plot(distances, E_CCSD_Np1_AB,
        label=r"cc-pVDZ cl-CCSD-eom, $\lambda$ = ?",
        color="darkred", linestyle='--',
        marker='^', linewidth=0.2,
        ms=2, mfc='none', mec='darkred', mew=0.3,
        markevery=2)

ax.plot(distances, E_CCSD_Nm1_AB,
        label=r"cc-pVDZ cl-CCSD-eom, $\lambda$ = ?",
        color="limegreen", linestyle='--',
        marker='v', linewidth=0.2,
        ms=2, mfc='none', mec='limegreen', mew=0.3,
        markevery=2)

ax.ticklabel_format(axis='y', style='plain', useOffset=False)

ax.legend(ncol=2, fontsize=4.5, frameon=False, loc='lower right')

plt.savefig("CCSD_v_Np1_v_Np2.pdf", dpi = 600)
