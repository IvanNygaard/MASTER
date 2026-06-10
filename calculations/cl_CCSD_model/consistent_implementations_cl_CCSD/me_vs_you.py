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


# Data from Alice and Ellen: 
alice_ellen_AB_aug_o_only = np.array([-152.530134027067, -152.578872834779, -152.612586124783, -152.636344103861, -152.653936641294,
                                       -152.666376649608, -152.675376955769, -152.681648359813, -152.686421437113, -152.689725263287,
                                       -152.692021734772, -152.693623772072, -152.694737493990, -152.695437509781, -152.695904335069,
                                       -152.696158220958, -152.696282709355, -152.696308876455, -152.696259061011, -152.696157192575,
                                       -152.696023296116, -152.695875873002, -152.695723711863, -152.695568995156, -152.695419489939,
                                       -152.695278054445, -152.695146698324, -152.695025932947, -152.694915643129, -152.694815087583,
                                       -152.694723077171, -152.694638887254, -152.694561512407, -152.694490136374, -152.694424073969,
                                       -152.694362677919, -152.694305617683, -152.694252375858, -152.694202726526, -152.694156374717,
                                       -152.694113103730, -152.693936771888, -152.693813040303, -152.693724883717, -152.693660250764,
                                       -152.693611681346, -152.693574444077, -152.693545391732, -152.693522374588, -152.693503892247,
                                       -152.693488873284])
 
alice_ellen_A_aug_o_only = np.array([-76.311156144570, -76.315530254255, -76.319626491733, -76.324831915338, -76.328607822082,
                                            -76.331708595254, -76.334315735261, -76.336618372344, -76.338553764540, -76.340172642259,
                                            -76.341532112449, -76.342597538374, -76.343456294108, -76.344163511777, -76.344717935227,
                                            -76.345158777805, -76.345510688566, -76.345786367875, -76.345997880381, -76.346153361731,
                                            -76.346268397884, -76.346356460698, -76.346424296968, -76.346473822627, -76.346511907919,
                                            -76.346542021993, -76.346567032119, -76.346588694257, -76.346607970933, -76.346625078585,
                                            -76.346640398505, -76.346653878957, -76.346665579997, -76.346675444912, -76.346683672664,
                                            -76.346690349769, -76.346695597247, -76.346699602802, -76.346702506724, -76.346704552087,
                                            -76.346705904497, -76.346707153469, -76.346706459995, -76.346706209045, -76.346706133265,
                                            -76.346706097205, -76.346706085196, -76.346706083369, -76.346706083147, -76.346706083262,
                                            -76.346706083912])
 
alice_ellen_B_aug_o_only = np.array([-76.337196066458, -76.338996239688, -76.340306688192, -76.340178139029, -76.342475924521,
                                            -76.342598089546, -76.343345071867, -76.343951278218, -76.344461299644, -76.344830109578,
                                            -76.345185820771, -76.345464342462, -76.345695663510, -76.345901594363, -76.346047599857,
                                            -76.346179467034, -76.346284617374, -76.346362039543, -76.346426824607, -76.346481493846,
                                            -76.346527815131, -76.346562055721, -76.346589856688, -76.346612045917, -76.346629610443,
                                            -76.346643472666, -76.346654341030, -76.346662803298, -76.346669349064, -76.346674371947,
                                            -76.346678185554, -76.346681025861, -76.346683193746, -76.346684868723, -76.346686115687,
                                            -76.346687038446, -76.346687771621, -76.346688329807, -76.346688751562, -76.346689067188,
                                            -76.346689303210, -76.346689810391, -76.346689871482, -76.346689918172, -76.346689941425,
                                            -76.346689950482, -76.346689951345, -76.346689953099, -76.346689953233, -76.346689953253,
                                            -76.346689951277])

alice_ellen_AB_aug_DZ = np.array([-152.386638972963, -152.432543118534, -152.466002427545, -152.489639818070, -152.506413748180,
                                       -152.518334072962, -152.526754920940, -152.532634688605, -152.536677709105, -152.539311031690,
                                       -152.541115991952, -152.542271664170, -152.542965130727, -152.543344631055, -152.543515136281,
                                       -152.543541519939, -152.543493823373, -152.543419793110, -152.543327391180, -152.543241806645,
                                       -152.543168628286, -152.543108203008, -152.543058615669, -152.543004911914, -152.542967397113,
                                       -152.542931789125, -152.542896218056, -152.542858881995, -152.542818795397, -152.542775281149,
                                       -152.542729457876, -152.542681040753, -152.542630609813, -152.542579191633, -152.542527450415,
                                       -152.542475881119, -152.542425452297, -152.542376671058, -152.542329878107, -152.542285614732,
                                       -152.542243788461, -152.542072233500, -152.541951002807, -152.541863067466, -152.541797915751,
                                       -152.541749011045, -152.541711824546, -152.541683065612,  -152.541660428715, -152.541642338180,
                                       -152.541627683152])
 
alice_ellen_A_aug_DZ = np.array([-76.246483203856, -76.249074097108, -76.250922615870, -76.254197862411, -76.257050749960,
                                            -76.259475263438, -76.261506596879, -76.263149090086, -76.264452288700, -76.265316140576,
                                            -76.266149839467, -76.266812568981, -76.267330160558, -76.267730202447, -76.268045854778,
                                            -76.268302633087, -76.268526877447, -76.268734222715, -76.268939265061, -76.269130632984,
                                            -76.269320952766, -76.269497722840, -76.269664567847, -76.269828474715, -76.269979032501,
                                            -76.270114459806, -76.270234729734, -76.270338721747, -76.270428825170, -76.270503041853,
                                            -76.270563760408, -76.270613524791, -76.270652743821, -76.270683288692, -76.270706514752,
                                            -76.270724015489, -76.270736816610, -76.270745992626, -76.270752603165, -76.270757070004,
                                            -76.270760127352, -76.270766035969, -76.270766786765, -76.270766018007, -76.270764938192,
                                            -76.270764137574, -76.270763692005, -76.270763467541, -76.270763358366, -76.270763304293,
                                            -76.270763279693])
 
alice_ellen_B_aug_DZ = np.array([-76.260432929400, -76.262819593810, -76.264296320650, -76.265426703086, -76.266332686864,
                                            -76.267073196945, -76.267709889914, -76.268248045818, -76.268733630400, -76.269209713483,
                                            -76.269554910897, -76.269844595036, -76.270090555942, -76.270287154116, -76.270448774380,
                                            -76.270571726438, -76.270661315620, -76.270723649571, -76.270765526706, -76.270788487169,
                                            -76.270808956153, -76.270817158571, -76.270823004685, -76.270825119690, -76.270825033794,
                                            -76.270823777986, -76.270821808708, -76.270819433054, -76.270816739495, -76.270814026377,
                                            -76.270811050035, -76.270807899208, -76.270804670780, -76.270801630895, -76.270798622781,
                                            -76.270795856750, -76.270793333302, -76.270791047144, -76.270789015240, -76.270787264388,
                                            -76.270785771911, -76.270781329860, -76.270779177265, -76.270777513571, -76.270776196746,
                                            -76.270775277873, -76.270774732856, -76.270774434782,  -76.270774283415, -76.270774212906,
                                            -76.270774182560])
 
 
E_INT_CL_CCSD_aug_cc_pVDZ_AE = (alice_ellen_AB_aug_DZ - alice_ellen_A_aug_DZ - alice_ellen_B_aug_DZ) *  627.5095                # kcal/mol
E_INT_CL_CCSD_aug_o_TZ_AE    = (alice_ellen_AB_aug_o_only - alice_ellen_A_aug_o_only - alice_ellen_B_aug_o_only) * 627.5095

  
# Plotting
fig, ax = plt.subplots()

plt.xlabel("Distance [Å]")
plt.ylabel(r"Interaction Energy [kcal mol$^{-1}$]")


ax.plot(distances, E_INT_CCSD_aug_cc_pVTZ, label = "aug-cc-pVTZ CCSD", color = "black", linestyle='-', linewidth=0.5)                             # reference 

ax.plot(distances, E_INT_CL_CCSD_aug_cc_pVDZ, label = "aug-cc-pVDZ cl-CCSD Impl. 1"               ,color = "limegreen", linestyle='-', linewidth=0.2, ms = 2, marker = 'd', mfc='none', mew=0.3, markevery=2)  
ax.plot(distances, E_INT_CL_CCSD_aug_cc_pVDZ_AE[10:], label = "aug-cc-pVDZ cl-CCSD Impl. 2"       , color = "darkblue", linestyle='-', linewidth=0.2, ms = 2, marker = '+',  mfc='none', mew=0.3, markevery=2)   
ax.plot(distances, E_INT_CL_CCSD_aug_o_TZ, label = "aug-cc-pVTZ / cc-pVTZ cl-CCSD Impl. 1"        , color = "red", linestyle='-', linewidth=0.2, ms = 2, marker = 'x',  mfc='none', mew=0.3, markevery=2)        
ax.plot(distances, E_INT_CL_CCSD_aug_o_TZ_AE[10:], label = "aug-cc-pVTZ / cc-pVTZ cl-CCSD Impl. 2", color = "purple", linestyle='-', linewidth=0.2, ms = 2, marker = 's',  mfc='none', mew=0.3, markevery=2)
           
ax.ticklabel_format(axis='y', style='plain', useOffset=False)


# Inset zoom (3–4 Å) showing only CL-CCSD curves
axins = inset_axes(ax, width="35%", height="35%", loc="center right")

axins.plot(distances, E_INT_CL_CCSD_aug_cc_pVDZ, label = "aug-cc-pVDZ cl-CCSD Impl. 1"               ,color = "limegreen", linestyle='-', linewidth=0.2, ms = 2, marker = 'd', mfc='none', mew=0.3, markevery=2)   
axins.plot(distances, E_INT_CL_CCSD_aug_cc_pVDZ_AE[10:], label = "aug-cc-pVDZ cl-CCSD Impl. 2"       , color = "darkblue", linestyle='-', linewidth=0.2, ms = 2, marker = '+',  mfc='none', mew=0.3, markevery=2)    
axins.plot(distances, E_INT_CL_CCSD_aug_o_TZ, label = "aug-cc-pVTZ / cc-pVTZ cl-CCSD Impl. 1"        , color = "red", linestyle='-', linewidth=0.2, ms = 2, marker = 'x',  mfc='none', mew=0.3, markevery=2)        
axins.plot(distances, E_INT_CL_CCSD_aug_o_TZ_AE[10:], label = "aug-cc-pVTZ / cc-pVTZ cl-CCSD Impl. 2", color = "purple", linestyle='-', linewidth=0.2, ms = 2, marker = 's',  mfc='none', mew=0.3, markevery=2)

axins.set_xlim(3.0, 4.5)
axins.set_ylim(-2.0, -0.5)  
axins.tick_params(labelsize=4)

mark_inset(ax, axins, loc1=1, loc2=3, fc="none", lw=0.4, color="black")

ax.legend(ncol=2, fontsize=3.5, frameon=False, loc='lower right')

plt.savefig("me_vs_them.pdf", dpi = 600)


# residual plot
fig, ax = plt.subplots()
plt.xlabel("q [Å]")
plt.ylabel(r"$\Delta E_{int}$ [kcal mol$^{-1}$]")
ax.scatter(distances, np.abs(E_INT_CL_CCSD_aug_o_TZ_AE[10:] - E_INT_CL_CCSD_aug_o_TZ), label = r"|$E_{int, Impl. 2} - E_{int, Impl. 1}$| aug-cc-pVTZ/cc-pVTZ", color = "black", s = 3, marker = "x")   
ax.legend(ncol=2, fontsize=5.5, frameon=False, loc='upper right')
plt.savefig("me_vs_them_residual_aug_TZ.pdf", dpi = 600)
ax.clear()
plt.xlabel("q [Å]")
plt.ylabel(r"$\Delta E_{int}$ [kcal mol$^{-1}$]")
ax.scatter(distances, np.abs(E_INT_CL_CCSD_aug_cc_pVDZ_AE[10:] - E_INT_CL_CCSD_aug_cc_pVDZ), label = r"|$E_{int, Impl. 2} - E_{int, Impl. 1}$| aug-cc-pVDZ", color = "black", s = 3, marker = "x")  
ax.legend(ncol=2, fontsize=5.5, frameon=False, loc='upper right')
plt.savefig("me_vs_them_residual_DZ.pdf", dpi = 600)

