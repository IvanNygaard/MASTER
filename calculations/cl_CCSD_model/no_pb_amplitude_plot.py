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
cc_pVDZ_data       = pd.read_csv("cc-pVDZ_CCSD_v_CL_CCSD.txt", sep=";")
cc_pVTZ_data       = pd.read_csv("cc-pVTZ_CCSD_v_CL_CCSD.txt", sep=";")
aug_cc_pVDZ_data   = pd.read_csv("aug-cc-pVDZ_CCSD_v_CL_CCSD.txt", sep=";")
aug_cc_pVTZ_data   = pd.read_csv("aug-cc-pVTZ_CCSD_v_CL_CCSD.txt", sep=";")
aug_o_only_TZ_data = pd.read_csv("aug-only-o-TZ.txt", sep=";")

# Remove whitespace
cc_pVDZ_data.columns       = cc_pVDZ_data.columns.str.strip()
cc_pVTZ_data.columns       = cc_pVTZ_data.columns.str.strip()
aug_cc_pVDZ_data.columns   = aug_cc_pVDZ_data.columns.str.strip()
aug_cc_pVTZ_data.columns   = aug_cc_pVTZ_data.columns.str.strip()
aug_o_only_TZ_data.columns = aug_o_only_TZ_data.columns.str.strip()


# Get results
# cc-pVDZ
distances   = cc_pVDZ_data["R"]                                                        # Å 

# Normal CCSD (with-cross)
E_CCSD         = cc_pVDZ_data["E_CCSD"]                                                # a.u. 
E_CCSD_A_AB    = cc_pVDZ_data["E_CCSD_A"]                                        
E_CCSD_B_AB    = cc_pVDZ_data["E_CCSD_B"]                                       
    
# Charge localized (no-cross)
E_CL_CCSD      = cc_pVDZ_data["E_CL_CCSD"]                                       
E_CL_CCSD_A_AB = cc_pVDZ_data["E_CL_CCSD_A"]                                    
E_CL_CCSD_B_AB = cc_pVDZ_data["E_CL_CCSD_B"]                                     

# Interaction energies
E_INT_CCSD_cc_pVDZ    = (E_CCSD    - E_CCSD_A_AB    - E_CCSD_B_AB)    * 627.5095       # kcal/mol 
E_INT_CL_CCSD_cc_pVDZ = (E_CL_CCSD - E_CL_CCSD_A_AB - E_CL_CCSD_B_AB) * 627.5095 


# cc-pVTZ
# Normal CCSD (with-cross)
E_CCSD         = cc_pVTZ_data["E_CCSD"]                                                # a.u.
E_CCSD_A_AB    = cc_pVTZ_data["E_CCSD_A"]                                       
E_CCSD_B_AB    = cc_pVTZ_data["E_CCSD_B"]                                        

# Charge localized (no-cross)
E_CL_CCSD      = cc_pVTZ_data["E_CL_CCSD"]                                      
E_CL_CCSD_A_AB = cc_pVTZ_data["E_CL_CCSD_A"]                                    
E_CL_CCSD_B_AB = cc_pVTZ_data["E_CL_CCSD_B"]                                    

# Interaction energies
E_INT_CCSD_cc_pVTZ    = (E_CCSD    - E_CCSD_A_AB    - E_CCSD_B_AB)    * 627.5095       # kcal/mol 
E_INT_CL_CCSD_cc_pVTZ = (E_CL_CCSD - E_CL_CCSD_A_AB - E_CL_CCSD_B_AB) * 627.5095 


# aug-cc-pVDZ
# Normal CCSD (with-cross)
E_CCSD         = aug_cc_pVDZ_data["E_CCSD"]                                            # a.u.
E_CCSD_A_AB    = aug_cc_pVDZ_data["E_CCSD_A"]    
E_CCSD_B_AB    = aug_cc_pVDZ_data["E_CCSD_B"]


# Charge localized (no-cross)
E_CL_CCSD      = aug_cc_pVDZ_data["E_CL_CCSD"]  
E_CL_CCSD_A_AB = aug_cc_pVDZ_data["E_CL_CCSD_A"] 
E_CL_CCSD_B_AB = aug_cc_pVDZ_data["E_CL_CCSD_B"] 


E_INT_CCSD_aug_cc_pVDZ    = (E_CCSD    - E_CCSD_A_AB    - E_CCSD_B_AB)     * 627.5095  # kcal/mol
E_INT_CL_CCSD_aug_cc_pVDZ = (E_CL_CCSD - E_CL_CCSD_A_AB - E_CL_CCSD_B_AB)  * 627.5095     


# aug-cc-pVTZ/cc-pVTZ
# Charge localized (no-cross)
E_CL_CCSD      = aug_o_only_TZ_data["E_CL_CCSD"]                                      # a.u.
E_CL_CCSD_A_AB = aug_o_only_TZ_data["E_CL_CCSD_A"] 
E_CL_CCSD_B_AB = aug_o_only_TZ_data["E_CL_CCSD_B"] 

# Interaction energies
E_INT_CL_CCSD_aug_o_TZ = (E_CL_CCSD - E_CL_CCSD_A_AB - E_CL_CCSD_B_AB) * 627.5095     # kcal/mol


# aug-cc-pVTZ
# Normal CCSD (with-cross)
E_CCSD         = aug_cc_pVTZ_data["E_CCSD"]                                            # a.u.
E_CCSD_A_AB    = aug_cc_pVTZ_data["E_CCSD_A"]   
E_CCSD_B_AB    = aug_cc_pVTZ_data["E_CCSD_B"]    

# Charge localized (no-cross)
E_CL_CCSD      = aug_cc_pVTZ_data["E_CL_CCSD"]   
E_CL_CCSD_A_AB = aug_cc_pVTZ_data["E_CL_CCSD_A"] 
E_CL_CCSD_B_AB = aug_cc_pVTZ_data["E_CL_CCSD_B"] 

# Interaction energies
E_INT_CCSD_aug_cc_pVTZ    = (E_CCSD    - E_CCSD_A_AB    - E_CCSD_B_AB)    * 627.5095  # kcal/mol
E_INT_CL_CCSD_aug_cc_pVTZ = (E_CL_CCSD - E_CL_CCSD_A_AB - E_CL_CCSD_B_AB) * 627.5095


# Plotting
fig, ax = plt.subplots()

plt.xlabel("Distance [Å]")
plt.ylabel(r"Interaction Energy [kcal mol$^{-1}$]")

ax.plot(distances, E_INT_CCSD_cc_pVDZ, label = "cc-pVDZ CCSD", color = "forestgreen", linestyle='-', linewidth=0.2, ms = 2, marker = 'd', mfc='none', mew=0.3, markevery=2)                                     # cc-pVDZ             CCSD
ax.plot(distances, E_INT_CCSD_aug_cc_pVDZ, label = "aug-cc-pVDZ CCSD", color = "darkblue", linestyle='-', linewidth=0.2, ms = 2, marker = '+',  mfc='none', mew=0.3, markevery=2)                             # aug-cc-pVDZ         CCSD
ax.plot(distances, E_INT_CCSD_aug_cc_pVTZ, label = "aug-cc-pVTZ CCSD", color = "red", linestyle='-', linewidth=0.2, ms = 2, marker = 'x',  mfc='none', mew=0.3, markevery=2)                               # aug-cc-pVTZ         CCSD
ax.plot(distances, E_INT_CL_CCSD_cc_pVDZ, label = "cc-pVDZ cl-CCSD", color = "forestgreen", linestyle='-', linewidth=0.2, ms = 2, marker = 's',  mfc='none', mew=0.3, markevery=2)                             # cc-pVDZ             CL-CCSD
ax.plot(distances, E_INT_CL_CCSD_aug_cc_pVDZ, label = "aug-cc-pVDZ cl-CCSD", color = "darkblue", linestyle='-', linewidth=0.2, ms = 2, marker = '^',  mfc='none', mew=0.3, markevery=2)                       # aug-cc-pVDZ         CL-CCSD
ax.plot(distances, E_INT_CL_CCSD_aug_o_TZ, label = "aug-cc-pVTZ / cc-pVTZ cl-CCSD", color = "red", linestyle = '-', linewidth=0.2, ms = 2, marker = 'v',  mfc='none', mew=0.3, markevery=2)                  # aug-cc-pVTZ/cc-pVTZ CL-CCSD

ax.ticklabel_format(axis='y', style='plain', useOffset=False)

# Inset zoom (3–4 Å) showing only CL-CCSD curves
axins = inset_axes(ax, width="35%", height="35%", loc="center right")

axins.plot(distances, E_INT_CL_CCSD_cc_pVDZ, color="forestgreen", linestyle='-', linewidth=0.2, ms = 2, marker = 's',  mfc='none', mew=0.3, markevery=2)
axins.plot(distances, E_INT_CL_CCSD_aug_cc_pVDZ, color="darkblue", linestyle='-', linewidth=0.2, ms = 2, marker = '^', mfc='none', mew=0.3, markevery=2)
axins.plot(distances, E_INT_CL_CCSD_aug_o_TZ, color="red", linestyle='-', linewidth=0.2, ms = 2, marker = 'v', mfc='none', mew=0.3, markevery=2)


axins.set_xlim(3.0, 4.5)
axins.set_ylim(-2.0, -0.5)  
axins.tick_params(labelsize=4)

mark_inset(ax, axins, loc1=1, loc2=3, fc="none", lw=0.4, color="black")

ax.legend(ncol=2, fontsize=4.5, frameon=False, loc='lower right')

plt.savefig("amplitudes_testing.pdf", dpi = 600)
