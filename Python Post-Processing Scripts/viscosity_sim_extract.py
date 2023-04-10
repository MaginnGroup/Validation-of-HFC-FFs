# -*- coding: utf-8 -*-

import numpy as np


# >>>>>>>>>>>>>>>>>> computed par atom temp profile <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

ext_final_profile = np.array([])

output_profile_final = open('velocity_profile.INDEX_viscosity.profile')
for line in output_profile_final:
  ext_profile = ((line.split(' ')[-1]))
  ext_final_profile = np.append(ext_final_profile, ext_profile)

ext_final_profile = np.array(ext_final_profile)
ext_final_profile = ext_final_profile[3:]
ext_final_profile = ext_final_profile.astype(float)
np.savetxt("ext_final_velocity_profile_INDEX_viscosity_profile.txt", ext_final_profile)