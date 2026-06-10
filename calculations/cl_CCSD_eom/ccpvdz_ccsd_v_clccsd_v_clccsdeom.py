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
cc_pVDZ_data       = pd.read_csv("cc-pVDZ_CCSD_v_CL_CCSD_v_CL_CCSD_EOM.txt", sep=";")


# Remove whitespace
cc_pVDZ_data.columns       = cc_pVDZ_data.columns.str.strip()


# Get results
# cc-pVDZ
distances   = cc_pVDZ_data["R"]                                                        # Å 


# Normal CCSD (with-cross)
E_CCSD         = cc_pVDZ_data["E_CCSD"]                                                # a.u. 
E_CCSD_A_AB    = cc_pVDZ_data["E_CCSD_A"]                                        
E_CCSD_B_AB    = cc_pVDZ_data["E_CCSD_B"]                                       


print(E_CCSD[0])
print(E_CCSD_A_AB[0])
print(E_CCSD_B_AB[0])



# Charge localized (no-cross)
E_CL_CCSD      = cc_pVDZ_data["E_CL_CCSD"]                                       
E_CL_CCSD_A_AB = cc_pVDZ_data["E_CL_CCSD_A"]                                    
E_CL_CCSD_B_AB = cc_pVDZ_data["E_CL_CCSD_B"]                                     


# Charge localized (no-cross)
E_CL_CCSD_EOM      = cc_pVDZ_data["E_CL_CCSD_EOM"]                                       
E_CL_CCSD_EOM_A_AB = cc_pVDZ_data["E_CL_CCSD_EOM_A"]                                    
E_CL_CCSD_EOM_B_AB = cc_pVDZ_data["E_CL_CCSD_EOM_B"]                                     


# Interaction energies
E_INT_CCSD_cc_pVDZ    = (E_CCSD    - E_CCSD_A_AB    - E_CCSD_B_AB)    * 627.5095       # kcal/mol 
E_INT_CL_CCSD_cc_pVDZ = (E_CL_CCSD - E_CL_CCSD_A_AB - E_CL_CCSD_B_AB) * 627.5095 
E_INT_CL_CCSD_EOM_cc_pVDZ = (E_CL_CCSD_EOM - E_CL_CCSD_EOM_A_AB - E_CL_CCSD_EOM_B_AB) * 627.5095 


# Plotting
fig, ax = plt.subplots()

plt.xlabel("Distance [Å]")
plt.ylabel(r"Interaction Energy [kcal mol$^{-1}$]")

ax.plot(distances, E_INT_CCSD_cc_pVDZ, label = "cc-pVDZ CCSD", color = "black", linestyle='-', linewidth=0.2, ms = 2, marker = 'x', mfc='none', mew=0.3, markevery=2)                     # cc-pVDZ             CCSD
ax.plot(distances, E_INT_CL_CCSD_cc_pVDZ, label = "cc-pVDZ cl-CCSD", color = "red", linestyle='-', linewidth=0.2, ms = 2, marker = 's',  mfc='none', mew=0.3, markevery=2)              # cc-pVDZ             CL-CCSD
ax.plot(distances, E_INT_CL_CCSD_EOM_cc_pVDZ, label = "cc-pVDZ cl-CCSD-eom", color = "forestgreen", linestyle = '-', linewidth=0.2, ms = 2, marker = 'd',  mfc='none', mew=0.3, markevery=2)            # cc-pVDZ             CL-CCSD-EOM


ax.ticklabel_format(axis='y', style='plain', useOffset=False)

ax.legend(ncol=2, fontsize=4.5, frameon=False, loc='lower right')

plt.savefig("first_eom_test.pdf", dpi = 600)
