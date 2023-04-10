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

def relative_abs_error(a,b):
  error = np.abs(a-b)
  rae = error/a
  return rae

def test_equilibration(running_average_data):  #5e-4
  start_index_frac = 0.2
  test_index_array = np.arange(start_index_frac, 0.8001, 0.005)
  test_index = test_index_array*len(running_average_data)
  test_points_index = test_index.astype(int)
  
  equilibration_point = 0
  for i in range(len(test_points_index)):
    check_index = test_points_index[i]
    equilibration_point = check_index
    for j in range(len(running_average_data[0,:])):
      check_range = np.arange(0, 2500.01, 1)
      check_range = check_range.astype(int)
      for k in check_range:
        if relative_abs_error(running_average_data[check_index,j], running_average_data[(check_index+k),j]) > 0.0005:
          break
      else:
        continue
      break
    else:
      break
  return equilibration_point, test_points_index

def swap_columns(df, col1, col2):
    col_list = list(df.columns)
    x, y = col_list.index(col1), col_list.index(col2)
    col_list[y], col_list[x] = col_list[x], col_list[y]
    df = df[col_list]
    return df

def running_average_array(array):
  running_ave = np.zeros([len(array)])
  for i in range(len(array)):
    mov_average = np.mean(array[0:i+1])
    running_ave[i] = mov_average
  return(running_ave)

ext_lo = np.loadtxt("ext_final_temp_profile_INDEX_thermalcond.txt")

ext_lo = np.array(ext_lo)

ext_lo = ext_lo.astype(float)

base_num_atoms = 1000

for i in range(len(ext_lo)):
  if ext_lo[i] > base_num_atoms:
    ext_lo[i] = 0

ext_lo = [s for s in ext_lo if s != 0]

lammps_input = pd.read_table("thermalsim.INDEX_thermalcond.txt", delimiter=" ")


x_length = lammps_input.iloc[0,-3]
z_length = lammps_input.iloc[0,-1]

dE_lo = 0.0025

samp_interval = 1000
samp_interval = int(samp_interval)

layers  = 20

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

ext_lo = np.reshape(ext_lo, (int(samp_sets), int(layers)))

ext_lo_table = pd.DataFrame(ext_lo, columns = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 
                                                   '10', '11', '12', '13', '14', '15', 
                                                   '16', '17', '18', '19', '20'])

z_interval = np.linspace(0.025, 0.975, 20)

z_spacings = z_interval * z_length 

z_slope = z_spacings[3:7]

def R_squared_value(x,y):
  model = LinearRegression()
  model.fit(x.reshape(-1,1), y)
  r_squared_lintest = model.score(x.reshape(-1,1), y)

  return r_squared_lintest

def tempvspos(pos, slope, intercept):
  return slope * pos + intercept


def slope_tempvsz(x, y):
  popt, pconv = curve_fit(tempvspos, x, y)
  slope = popt[0]
  return slope


def kappa(inv_temp_grad, Q):
  kappa_real = -Q * inv_temp_grad
  kappa_realtosi = 69476948.46        # Conversion factor determined via units conversion analysis
  kappa_ = kappa_real * kappa_realtosi
  return kappa_

def block_sizes(x):
   block_size_list = np.array([])
   max_block_size = int(x/2)
   for i in range(1, (max_block_size + 1)):
     if max_block_size % i == 0:
       block_size_list = np.append(block_size_list, i)
   block_size_list = block_size_list.astype(int)
   return block_size_list

def kappa_ave_error(file1, steady_state_time, total_time_1, tblock, tgap, sampling_interval, dE):
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

  ave_temps = np.zeros([len(start_block_timestep_index), 10]) # 10 is total number of layers from hot to cold
  for i in range(len(start_block_timestep_index)):
    ave_temps[i,:] = file1.iloc[start_block_timestep_index[i]:(stop_block_timestep_index[i] + 1),0:10].mean(axis =0) #*******************
  ave_temps_slope = ave_temps[:,3:7]                                                                              
  dT_dz = np.zeros([len(start_block_timestep_index)])
  for v in range(len(start_block_timestep_index)):
    dT_dz[v] = slope_tempvsz(z_slope, ave_temps_slope[v,:])
  dz_dT = np.zeros([len(dT_dz)])
  for u in range(len(dT_dz)):
    dz_dT[u] = 1/dT_dz[u]
  R_squared_kappa = np.zeros([len(dT_dz)])
  for r in range(len(dT_dz)):
    R_squared_kappa[r] = R_squared_value(z_slope, ave_temps_slope[r,:])
  for w in range(len(R_squared_kappa)):
    if R_squared_kappa[w] < 0.95:
      dz_dT[w] = 0
    else:
      dz_dT[w] = dz_dT[w]
  dz_dT = [s for s in dz_dT if s != 0]
  dz_dT = np.array(dz_dT)
  dz_dT = dz_dT.astype(float)
  dQ_file1 = (dE * total_time_1)/(total_time_1 * x_length**2)
  kappa_file1 = kappa(dz_dT, dQ_file1)
  kappa_file1_ave = np.mean(kappa_file1)
  kappa_file1_var_unnormalized = np.var(kappa_file1)
  kappa_file1_var_normalized = kappa_file1_var_unnormalized/(num_blocks - 1)
  kappa_file1_std = np.sqrt(kappa_file1_var_normalized)
  kappa_file1_percentstd = kappa_file1_std/kappa_file1_ave * 100
  return(kappa_file1, kappa_file1_ave, kappa_file1_std, kappa_file1_percentstd, R_squared_kappa, dQ_file1, num_blocks)

def error_vs_blocks_data(block_sizes):
  tblocks = block_sizes * samp_interval
  #tblocks = np.array(tblocks)
  #tblocks = tblocks.reshape(-1,1)
  block_percent_error = np.array([])
  block_abs_error = np.array([])
  num_blocks  = np.array([])
  kappa_values = np.array([])
  max_number_blocks = (total_time - steady_state_time)/(np.min(tblocks))
  max_number_blocks = int(max_number_blocks)
  R_squared_values_array = np.zeros([len(block_sizes), max_number_blocks])
 
  for i in range(len(tblocks)):
    kappa_array = kappa_ave_error(ext_lo_table, steady_state_time, total_time, tblocks[i], 0, samp_interval, dE_lo)
    block_percent_error = np.append(block_percent_error, kappa_array[3])
    block_abs_error = np.append(block_abs_error, kappa_array[2])
    num_blocks =  np.append(num_blocks, kappa_array[-1])
    kappa_values =  np.append(kappa_values, kappa_array[1])
    R_squared_values_array[i,0:len(kappa_array[-3])] = kappa_array[-3]
  
  max_error_index_array = np.where(block_percent_error == np.amax(block_percent_error))
  max_error_index = max_error_index_array[0]
  max_error_index = int(max_error_index)
  final_kappa = kappa_values[max_error_index]
  final_abs_error = block_abs_error[max_error_index]
  final_percent_error = block_percent_error[max_error_index]
  R_squared_final = R_squared_values_array[max_error_index, :]
  R_squared_final = [s for s in R_squared_final if s != 0]
  final_tblock = tblocks[max_error_index]

  return(tblocks, num_blocks, block_percent_error, block_abs_error, kappa_values, final_tblock, R_squared_final, final_kappa, final_abs_error, final_percent_error) 
  #return(tblocks, num_blocks, block_percent_error, block_abs_error, kappa_values)



blocks = block_sizes(int(samp_sets_prod))

error_analysis = error_vs_blocks_data(blocks)

kappa_value = error_analysis[-3]
kappa_value = np.array([kappa_value])

kappa_value_abserror = error_analysis[-2]
kappa_value_abserror = np.array([kappa_value_abserror])

kappa_value_percenterror = error_analysis[-1]
kappa_value_percenterror = np.array([kappa_value_percenterror])

np.savetxt("kappa_value_INDEX.txt", kappa_value)
np.savetxt("kappa_value_percenterror_INDEX.txt", kappa_value_percenterror)
np.savetxt("kappa_value_abserror_INDEX.txt", kappa_value_abserror)
