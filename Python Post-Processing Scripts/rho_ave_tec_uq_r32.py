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
  rae = error/a
  return rae

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

Bridgman_Output_data = pd.read_table("Output_INDEX1.txt", sep='\s+')

Bridgman_Output_data = swap_columns(Bridgman_Output_data, 'P', 'V')

equilibration_start = np.max(np.loadtxt("tequib.txt"))
equilibration_start = int(equilibration_start)

final_equib_time = equilibration_start

total_time = int(15000000)

def time_blocksize_rho(total_time, final_equib_time, num_rho_sets):
  post_equib_time = total_time - final_equib_time
  if post_equib_time%num_rho_sets == 0:
    time_block_rho = post_equib_time/num_rho_sets
  elif (post_equib_time-1)%num_rho_sets == 0:
    time_block_rho = (post_equib_time-1)/num_rho_sets
  time_block_rho = int(time_block_rho)
  return time_block_rho

num_rho_sets = 3
time_blocks_rho = time_blocksize_rho(total_time, final_equib_time, num_rho_sets)

ave_rho_array = np.zeros(num_rho_sets)

for i in range(num_rho_sets):
  ave_rho = running_average_array(Bridgman_Output_data.iloc[int((final_equib_time+time_blocks_rho*i)/1000):int((final_equib_time+time_blocks_rho*(i+1))/1000),2])
  ave_rho_array[i] = ave_rho[-1]

np.savetxt('ave_density_array_INDEX1.txt', ave_rho_array)