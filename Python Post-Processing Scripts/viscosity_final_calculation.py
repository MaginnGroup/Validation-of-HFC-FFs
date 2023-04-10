# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd





from scipy.optimize import curve_fit

from sklearn.linear_model import LinearRegression


def running_average(file):
  running_ave = np.zeros([len(file.iloc[:,0]), len(file.iloc[0,:])])
  for i in range(len(file.iloc[0,:])):
    for j in range(len(file.iloc[:,0])):
      mov_average = file.iloc[0:j+1,i].mean() 
      running_ave[j,i] = mov_average
  return(running_ave)

visco_data = pd.read_table("momentumflux.INDEX_viscosity.txt", delimiter=" ")

total_momentum_flux = visco_data.iloc[-1,-1]

#print(total_momentum_flux)

ext_lo = np.loadtxt("ext_final_velocity_profile_INDEX_viscosity_profile.txt")

ext_lo = np.array(ext_lo)

ext_lo = ext_lo.astype(float)

base_num_atoms = 1000

for i in range(len(ext_lo)):
  if ext_lo[i] > base_num_atoms:
    ext_lo[i] = 0

ext_lo = [s for s in ext_lo if s != 0]

#print(len(ext_lo))



x_length = visco_data.iloc[-1,-4]

#print(x_length)

samp_interval = 1000
samp_interval = int(samp_interval)

layers  = 20

aspect_ratio = 20/3 

z_length = aspect_ratio * x_length

ext_start_time = 0

total_time = ((len(ext_lo))/layers) * samp_interval

total_time = int(total_time)

samp_sets = int(total_time/samp_interval)

samp_sets_edit = []

if samp_sets % 100 == 0:
  steady_state_time = int(total_time/2)
else: 
  samp_sets_edit = samp_sets + 100 - (samp_sets % 100)
  steady_state_time = total_time - ((samp_sets_edit * samp_interval)/2)


samp_sets_prod = int((total_time - steady_state_time)/samp_interval)

#print(samp_sets)

#print(samp_sets_prod)

ext_lo = np.reshape(ext_lo, (int(samp_sets), int(layers)))

ext_lo_table = pd.DataFrame(ext_lo, columns = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 
                                                   '10', '11', '12', '13', '14', '15', 
                                                   '16', '17', '18', '19', '20'])

time_steps = np.linspace(ext_start_time + samp_interval, total_time, int(samp_sets))



z_interval = np.linspace(0.025, 0.975, 20)

z_spacings = z_interval * z_length 

z_slope = z_spacings[11:20]

'''
def running_average(file):
  running_ave = np.zeros(len(file))
  for i in range(len(file)):
    mov_average = np.mean(file[0:i+1])
    running_ave[i] = mov_average
  return(running_ave)
'''

def R_squared_value(x,y):
  model = LinearRegression()
  model.fit(x.reshape(-1,1), y)
  r_squared_lintest = model.score(x.reshape(-1,1), y)

  return r_squared_lintest

def velvspos(pos, slope, intercept):
  return slope * pos + intercept

def slope_velvsz(x, y):
  popt, pconv = curve_fit(velvspos, x, y)
  slope = popt[0]
  return slope

def nu(inv_velocity_grad, J):
  nu_real = np.abs(J * inv_velocity_grad)
  nu_realtosi = 1.660539067e-2         # Conversion factor determined via units conversion analysis
  nu_ = nu_real * nu_realtosi
  return nu_

def block_sizes(x):
   block_size_list = np.array([])
   max_block_size = int(x/2)
   for i in range(1, (max_block_size + 1)):
     if max_block_size % i == 0:
       block_size_list = np.append(block_size_list, i)
   block_size_list = block_size_list.astype(int)
   return block_size_list

def nu_ave_error(file1, steady_state_time, total_time_1, tblock, tgap, sampling_interval, total_mom_flux):
  num_ave_points = tblock/sampling_interval
  num_blocks = (total_time_1 - steady_state_time)/(tblock + tgap)
  if num_blocks - int(num_blocks) > 0.5:
    num_blocks = num_blocks + 1
  else:
    num_blocks = num_blocks
  num_blocks = int(num_blocks)
  inv_samp_interval = 1/sampling_interval
  start_block_times = np.linspace(steady_state_time, total_time_1 - (tblock + tgap), num_blocks)
  start_block_timestep_index = start_block_times * inv_samp_interval
  start_block_timestep_index = start_block_timestep_index.astype(int)
  stop_block_times = start_block_times + tblock
  stop_block_timestep_index = stop_block_times * inv_samp_interval
  stop_block_timestep_index = stop_block_timestep_index.astype(int)


  ave_vel = np.zeros([len(start_block_timestep_index), 9])
  for i in range(len(start_block_timestep_index)):
    ave_vel[i,:] = file1.iloc[start_block_timestep_index[i]:(stop_block_timestep_index[i] + 1),11:20].mean(axis =0) #*******************
  ave_vel_slope = ave_vel[:,0:9]                                                                              
  dV_dz = np.zeros([len(start_block_timestep_index)])
  for v in range(len(start_block_timestep_index)):
    dV_dz[v] = slope_velvsz(z_slope, ave_vel_slope[v,:])
  dz_dV = np.zeros([len(dV_dz)])
  for u in range(len(dV_dz)):
    dz_dV[u] = 1/dV_dz[u]
  R_squared_nu = np.zeros([len(dV_dz)])
  for r in range(len(dV_dz)):
    R_squared_nu[r] = R_squared_value(z_slope, ave_vel_slope[r,:])
  for w in range(len(R_squared_nu)):
    if R_squared_nu[w] < 0.95:
      dz_dV[w] = 0
    else:
      dz_dV[w] = dz_dV[w]
  dz_dV = [s for s in dz_dV if s != 0]
  dz_dV = np.array(dz_dV)
  dz_dV = dz_dV.astype(float)
  dJ_file1 = (np.abs(total_mom_flux))/(2 * total_time_1 * x_length**2)
  nu_file1 = nu(dz_dV, dJ_file1)
  nu_file1_ave = np.mean(nu_file1)
  nu_file1_var_unnormalized = np.var(nu_file1)
  nu_file1_var_normalized = nu_file1_var_unnormalized/(num_blocks - 1)
  nu_file1_std = np.sqrt(nu_file1_var_normalized)
  nu_file1_percentstd = nu_file1_std/nu_file1_ave * 100
  return(nu_file1, nu_file1_ave, nu_file1_std, nu_file1_percentstd, R_squared_nu, dJ_file1, num_blocks)

def error_vs_blocks_data(block_sizes):
  tblocks = block_sizes * samp_interval
  #tblocks = np.array(tblocks)
  #tblocks = tblocks.reshape(-1,1)
  block_percent_error = np.array([])
  block_abs_error = np.array([])
  num_blocks  = np.array([])
  nu_values = np.array([])
  max_number_blocks = (total_time - steady_state_time)/(np.min(tblocks))
  max_number_blocks = int(max_number_blocks)
  R_squared_values_array = np.zeros([len(block_sizes), max_number_blocks])
 
  for i in range(len(tblocks)):
    nu_array = nu_ave_error(ext_lo_table, steady_state_time, total_time, tblocks[i], 0, samp_interval, total_momentum_flux)
    block_percent_error = np.append(block_percent_error, nu_array[3])
    block_abs_error = np.append(block_abs_error, nu_array[2])
    num_blocks =  np.append(num_blocks, nu_array[-1])
    nu_values =  np.append(nu_values, nu_array[1])
    R_squared_values_array[i,0:len(nu_array[-3])] = nu_array[-3]
  
  max_error_index_array = np.where(block_percent_error == np.amax(block_percent_error))
  max_error_index = max_error_index_array[0]
  max_error_index = int(max_error_index)
  final_nu = nu_values[max_error_index]
  final_abs_error = block_abs_error[max_error_index]
  final_percent_error = block_percent_error[max_error_index]
  R_squared_final = R_squared_values_array[max_error_index, :]
  R_squared_final = [s for s in R_squared_final if s != 0]
  final_tblock = tblocks[max_error_index]

  return(tblocks, num_blocks, block_percent_error, block_abs_error, nu_values, final_tblock, R_squared_final, final_nu, final_abs_error, final_percent_error) 
  #return(tblocks, num_blocks, block_percent_error, block_abs_error, nu_values)

blocks = block_sizes(int(samp_sets_prod))

error_analysis = error_vs_blocks_data(blocks)

#print('[',error_analysis[-3], ',' ,error_analysis[-2], ',',error_analysis[-1],']')

viscosity_value = error_analysis[-3]
viscosity_value = np.array([viscosity_value])

viscosity_value_abserror = error_analysis[-2]
viscosity_value_abserror = np.array([viscosity_value_abserror])

viscosity_value_percenterror = error_analysis[-1]
viscosity_value_percenterror = np.array([viscosity_value_percenterror])

np.savetxt("viscosity_value_INDEX.txt", viscosity_value)
np.savetxt("viscosity_value_percenterror_INDEX.txt", viscosity_value_percenterror)
np.savetxt("viscosity_value_abserror_INDEX.txt", viscosity_value_abserror)
