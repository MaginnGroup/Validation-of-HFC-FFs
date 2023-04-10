# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

def running_average(file):
  running_ave = np.zeros([len(file.iloc[:,0]), len(file.iloc[0,:])])
  for i in range(len(file.iloc[0,:])):
    for j in range(len(file.iloc[:,0])):
      mov_average = file.iloc[0:j+1,i].mean() 
      running_ave[j,i] = mov_average
  return(running_ave)

def relative_abs_error(a,b):
  error = np.abs(a-b)
  rae = np.abs(error/a)
  return rae

def test_equilibration(running_average_data, run_ave_error_tol, test_time_steps, start_index_frac, stop_index_frac):  #5e-4
  test_index_array = np.arange(start_index_frac, stop_index_frac+0.0001, 0.005)
  test_index = test_index_array*len(running_average_data)
  test_points_index = test_index.astype(int)
  equilibration_point = 0
  for i in range(len(test_points_index)):
    check_index = test_points_index[i]
    equilibration_point = check_index
    for j in range(len(running_average_data[0,:])):
      check_range = np.arange(0, test_time_steps+0.01, 1)
      check_range = check_range.astype(int)
      for k in check_range:
        if relative_abs_error(running_average_data[check_index,j], running_average_data[(check_index+k),j]) > run_ave_error_tol: # 0.0005:
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

Bridgman_Output_data = pd.read_table("Output_INDEX1_NVT_constrho.txt", sep='\s+')

Bridgman_Output_data = swap_columns(Bridgman_Output_data, 'P', 'V')

run_ave_Bridgman_data = running_average(Bridgman_Output_data.iloc[:,:])

equilibration_start = test_equilibration(run_ave_Bridgman_data[:,2:], 0.0005, 2500, 0.2, 0.8)

def get_final_equib_time(equib_point):
  equilibration_time_index = equib_point     
  rem = equilibration_time_index%500
  to_add = 0
  if rem != 0:
    to_add = 500 - rem
  final_equib_time = equilibration_time_index + to_add
  final_equib_time_ = int(final_equib_time)
  return final_equib_time_

final_equib_time = get_final_equib_time(equilibration_start[0])

final_equib_time_stored = np.array([final_equib_time])*1000
np.savetxt('final_equib_time_INDEX1.txt', final_equib_time_stored)

Bridgman_equib_data = running_average(Bridgman_Output_data.iloc[final_equib_time:,:])

ave_density = Bridgman_equib_data[-1,2]
ave_density = np.array([ave_density])
np.savetxt('ave_density_INDEX1.txt', ave_density)

kb = 0.0019872067        # kB in Kcal/mole/K
kb_si = 1.380649e-23  # kb in si units
R   = 8.3144626733792 # J/mol/K
kcal_per_mole_2_Joules = 6.947695411e-21   # using 1 kcal = 4184 J # 1 kcal/mol = 4184/Avogadr_num J
kcal_2_J = 4184
Angstrom_cube_2_m_cube = 1e-30
Avg_num = 6.0221408e23
atm2Pa = 101325
bar2Pa = 100000

N = 2000
Molar_Energy_real2SI = kcal_2_J/N

U_config = (Bridgman_Output_data.iloc[:,11] + Bridgman_Output_data.iloc[:,16] + Bridgman_Output_data.iloc[:,17]) * Molar_Energy_real2SI # Molar energy in SI
Pressure_Vol = (Bridgman_Output_data.iloc[:,1] * Bridgman_Output_data.iloc[:,4]) * \
                    Angstrom_cube_2_m_cube * (Avg_num/N) * atm2Pa
H_config = U_config + Pressure_Vol
ave_Hconfig = running_average_array(H_config)

Bridgman_Output_data['U_config'] = U_config

Bridgman_Output_data['H_config'] = H_config

del Bridgman_Output_data["Evdwl_intra"]

#print(Bridgman_Output_data.head())

def block_sizes(x):
   block_size_list = np.array([])
   max_block_size = int(x/2)
   for i in range(1, (max_block_size + 1)):
     if max_block_size % i == 0:
       block_size_list = np.append(block_size_list, i)
   block_size_list = block_size_list.astype(int)
   return block_size_list

def Bridgman_data_UQ(file1, equib_time, total_time_1, tblock, tgap, sampling_interval):
  num_ave_points = tblock/sampling_interval
  num_blocks = (total_time_1 - equib_time)/(tblock + tgap)
  if num_blocks - int(num_blocks) > 0.5:
    num_blocks = num_blocks + 1
  else:
    num_blocks = num_blocks
  num_blocks = int(num_blocks)
  inv_samp_interval = 1/sampling_interval
  start_block_times = np.linspace(equib_time, total_time_1 - (tblock + tgap), num_blocks)
  start_block_timestep_index = start_block_times * inv_samp_interval
  start_block_timestep_index = start_block_timestep_index.astype(int)
  stop_block_times = start_block_times + tblock
  stop_block_timestep_index = stop_block_times * inv_samp_interval
  stop_block_timestep_index = stop_block_timestep_index.astype(int)
  ave_data = np.zeros([len(start_block_timestep_index), len(file1.iloc[0,:])]) 
  for i in range(len(start_block_timestep_index)):
    ave_data[i,:] = file1.iloc[start_block_timestep_index[i]:(stop_block_timestep_index[i] + 1),:].mean(axis =0) #*******************
  
  data_file1_ave = np.mean(ave_data, axis=0)
  data_file1_var_unnormalized = np.var(ave_data, axis=0)
  data_file1_var_normalized = data_file1_var_unnormalized/(num_blocks - 1)
  data_file1_std = np.sqrt(data_file1_var_normalized)
  data_file1_percentstd = data_file1_std/(np.abs(data_file1_ave)) * 100
  return(ave_data, data_file1_ave, data_file1_std, data_file1_percentstd, num_blocks)



def error_vs_blocks_data(file1, block_sizes, sampling_interval, total_time, equib_time):
  tblocks = block_sizes * sampling_interval
  block_percent_error = np.zeros([len(block_sizes), len(file1.iloc[0,:])])
  block_abs_error = np.zeros([len(block_sizes), len(file1.iloc[0,:])])
  num_blocks  = np.zeros([len(block_sizes)])
  data_values = np.zeros([len(block_sizes), len(file1.iloc[0,:])])
  max_number_blocks = (total_time - equib_time)/(np.min(tblocks))
  max_number_blocks = int(max_number_blocks)
  for i in range(len(tblocks)):
    data_array = Bridgman_data_UQ(file1, equib_time, total_time, tblocks[i], 0, sampling_interval)
    block_percent_error[i,:] = data_array[3]
    block_abs_error[i,:] = data_array[2]
    num_blocks[i] =  data_array[-1]
    data_values[i,:] =  data_array[1]
  final_Bridgman_values = np.zeros([len(file1.iloc[0,:])])
  final_abs_error = np.zeros([len(file1.iloc[0,:])]) 
  final_percent_error = np.zeros([len(file1.iloc[0,:])]) 
  final_numblock = np.zeros([len(file1.iloc[0,:])]) 
  for j in range(len(file1.iloc[0,:])):      
    max_error_index_array = np.where(block_percent_error[:,j] == np.amax(block_percent_error[:,j]))
    max_error_index = max_error_index_array[0]
    #print(max_error_index)
    max_error_index = int(max_error_index)
    final_Bridgman_values[j] = data_values[max_error_index,j]
    final_abs_error[j] = block_abs_error[max_error_index,j]
    final_percent_error[j] = block_percent_error[max_error_index,j]
    final_numblock[j] = num_blocks[max_error_index]
  
  return(tblocks, num_blocks, block_percent_error, block_abs_error, data_values, final_numblock, final_Bridgman_values, final_abs_error, final_percent_error)

sampling_interval = 1000
total_time = len(Bridgman_Output_data.iloc[:,0]) * sampling_interval 
equilibrated_time = final_equib_time * sampling_interval 
samp_sets_prod = int((total_time - equilibrated_time)/sampling_interval)
blocks = block_sizes(samp_sets_prod)

error_analysis = error_vs_blocks_data(Bridgman_Output_data.iloc[:,2:], blocks, sampling_interval, total_time, equilibrated_time)

np.savetxt('Bridgman_ave_values_INDEX1.txt', error_analysis[-3][0:])
np.savetxt('Bridgman_abs_error_INDEX1.txt', error_analysis[-2][0:])
np.savetxt('Bridgman_percent_error_INDEX1.txt', error_analysis[-1][0:])
np.savetxt('Bridgman_num_blocks_INDEX1.txt', error_analysis[-4][0:])

np.savetxt('ave_Hconfig_INDEX1.txt', np.array([error_analysis[-3][-1]]))
np.savetxt('ave_Hconfig_percent_error_INDEX1.txt', np.array([error_analysis[-1][-1]]))

np.savetxt('ave_Uconfig_INDEX1.txt', np.array([error_analysis[-3][-2]]))
np.savetxt('ave_Uconfig_percent_error_INDEX1.txt', np.array([error_analysis[-1][-2]]))

np.savetxt('rho_INDEX1.txt', np.array([error_analysis[-3][0]]))
np.savetxt('rho_percent_error_INDEX1.txt', np.array([error_analysis[-1][0]]))



equilibration_start_pressure = test_equilibration(run_ave_Bridgman_data[:,1:3], 0.01, 4999, 0.2, 0.75)
#print(equilibration_start_pressure[0])

final_equib_time_press = get_final_equib_time(equilibration_start_pressure[0])
sampling_interval = 1000
total_time = len(Bridgman_Output_data.iloc[:,0]) * sampling_interval 
equilibrated_time_press = final_equib_time_press * sampling_interval 
samp_sets_prod_press = int((total_time - equilibrated_time_press)/sampling_interval)
blocks_press = block_sizes(samp_sets_prod_press)

final_equib_time_stored_press = np.array([final_equib_time_press])*1000
np.savetxt('final_equib_time_press_INDEX1.txt', final_equib_time_stored_press)

error_analysis_pressure = error_vs_blocks_data(Bridgman_Output_data.iloc[:,1:3], blocks_press, sampling_interval, total_time, equilibrated_time_press)

np.savetxt('ave_Press_INDEX1.txt', np.array([error_analysis_pressure[-3][0]]))
np.savetxt('ave_Press_percent_error_INDEX1.txt', np.array([error_analysis_pressure[-1][0]]))

