
This repository provides one or more example(s) of all relevant
LAMMPS input scripts,
Data files,
LAMMPS output files, and
Python post-processing scripts.

In these example scripts, all occurences of 'INDEX', 'INDEX1', 'DIR_INDEX' 
or 'DATA_INDEX' are simply placeholders to be substituted with relevant names to give actual file names.
In the same vein, 'TEMP_INDEX' and 'PRESS_INDEX' in the example LAMMPS input scripts are placeholders 
to be replaced with the appropriate values of temperature and pressure respectively.

Also, all mention or references to v52, v37, v36 and v23 for R32 are to be understood as versions a, b, c and d 
respectively as used in the paper for this work. v100, v80, v74 and v7 represent versions a, b, c and d 
respectively for R125. 

All files ending with the .data extension are simply FF data input files for simulations in LAMMPS. 

Using these example scripts for the large number of simulations performed requires some automation scripts 
which can be written using Bash scripting. Those Bash scripts used for job submission automation are not included here.
The files 'Labels_R32_NPTBridgman_allpoints.txt' and 'Labels_R32_Bridgman_NVT_allpoints.txt' 
contain the directory (and file) labels used for the automations. 

One of those labels is 'R32_v52_2000mol_nptBridgman_273K_273P_realstate'. If all occurences of 
'INDEX', 'INDEX1', and/or 'DIR_INDEX' are replaced with 'R32_v52_2000mol_nptBridgman_273K_273P_realstate', and 'DATA_INDEX' replaced with R32_v52_2000mol, then the example NPT and NVT LAMMPS input scripts can be run to give the output files; Output_R32_v52_2000mol_nptBridgman_273K_273P_realstate.txt 
and Output_R32_v52_2000mol_nptBridgman_273K_273P_realstate_NVT_Cv_statepoint.txt
which have been provided as example LAMMPS output files from NPT and NVT simulations respectively.


Likewise,
momentumflux.R32_v52_2000mol_nptBridgman_273K_273P_realstate_viscosity.txt,
velocity_profile.R32_v52_2000mol_nptBridgman_273K_273P_realstate_viscosity.profile,
temp_profile.R32_v52_2000mol_nptBridgman_273K_273P_realstate_thermalcond.0.0025,
massdens_profile.R32_v52_2000mol_nptBridgman_273K_273P_realstate_thermalcond.0.0025,
lammps_temp_profile.R32_v52_2000mol_nptBridgman_273K_273P_realstate_thermalcond.0.0025,
numdens_profile.R32_v52_2000mol_nptBridgman_273K_273P_realstate_thermalcond.0.0025, and
thermalsim.R32_v52_2000mol_nptBridgman_273K_273P_realstate_thermalcond
are example LAMMPS output files obtained from running the in.lmp.INDEX_reshape3r20.txt followed
by the in.lmp.INDEX_thermalcond.txt and in.lmp.INDEX_viscosity.txt LAMMPS input scripts with 
the label placeholders replaced with 'R32_v52_2000mol_nptBridgman_273K_273P_realstate'.


The Python post-processing scripts
bridgman_calculations_uncertainty_quantification_phase1.py , 
rho_ave_tec_uq_r32.py , and
bridgman_calculations_uq_nvt_phase1.py
Take in the appropriate LAMMPS output files and give average values of relevant quantities (with uncertainties) 
required for final calculations of all thermodynamic properties from NPT and NVT simulations. 
These quantities include the configurational enthalpies, internal energies, densities and pressure.

These scripts can be tested using the NPT and NVT output files provided for the simulation
with label 'R32_v52_2000mol_nptBridgman_273K_273P_realstate' .

Psi4_input_R32_273K.txt  and Psi4_R32_273K_output.txt are example input and output files for
the QM calculations for R32. Those for R125 are also available.

The Python post-processing scripts
thermal_sim_extract.py, and viscosity_sim_extract.py are used for preliminary 
post-processing of the LAMMPS output files from thermal conductivity and viscosity simulations respectively.
Example output files from these preliminary post-processing for 'R32_v52_2000mol_nptBridgman_273K_273P_realstate' are;
ext_final_number_density_profile_R32_v52_2000mol_nptBridgman_273K_273P_realstate_thermalcond.txt,
ext_final_mass_density_profile_R32_v52_2000mol_nptBridgman_273K_273P_realstate_thermalcond.txt,
ext_final_lammps_temp_profile_R32_v52_2000mol_nptBridgman_273K_273P_realstate_thermalcond.txt,
ext_final_temp_profile_R32_v52_2000mol_nptBridgman_273K_273P_realstate_thermalcond.txt, and
ext_final_velocity_profile_R32_v52_2000mol_nptBridgman_273K_273P_realstate_viscosity_profile.txt

Then the Python scripts
thermal_cond_finalcalc.py, and viscosity_final_calculation.py 
can be used to further process these preliminary post-processing outputs to get final average values (with uncertainties) of 
thermal conductivty and viscosity for a single state point and FF version.

The outlined process using 'R32_v52_2000mol_nptBridgman_273K_273P_realstate' as an example 
was completed for all other simulation points as detailed in the Label files; 'Labels_R32_NPTBridgman_allpoints.txt' 
and 'Labels_R32_Bridgman_NVT_allpoints.txt'.

Then, the Jupyter notebooks R32_Bridgman_final_calcs_et_plots_ave_rho_array_TEC.ipynb , and
R125_Bridgman_final_calcs_et_plots_ave_rho_array_TEC.ipynb can be used to perform final calculations
of the thermodynamic properties and generate the plots and result tables in the paper. 
These notebooks also take in the already completely calculated thermal conductivities and
viscosities and generate the relevant plots.

Calculations of self-diffusivity and RDFs use output NVT trajectory files and the PyLAT package. More information on using PyLAT along with the necessary scripts and codes are available at https://github.com/MaginnGroup/PyLAT .



