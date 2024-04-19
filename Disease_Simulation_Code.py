#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 08:57:04 2021

@author: melissacollier

This code will run a percolation disease simulation over a range of T values for all networks in list "Synchrony Networks"
As downloaded, this code will run disease simulations on 25 networks modeled from PC data on three T values
"""
import os

import numpy as np

import Simulation_Functions as sim
import pandas as pd


Synchrony_Networks =sim.get_networks_with_directory_specifc('generated_graphs', ".graphml")

num_sims = 100


synch_size = []
synch_prob = []
synch_AM = []
synch_JX = []
synch_AF= []
AM_total = []
AF_total = []
JX_total = []
R0s_list = []
graph_list = []
T_list_f = []
AMAX_list = []
AFAX_list = []


total_graphs = len(Synchrony_Networks)

T_list = [0.03, 0.055, 0.075] #Transmissibility values to try for the disease simulation
i = 1
for G in Synchrony_Networks:
    for t in T_list:
        print("Working on T ",t, " sims for synch graph ", i, " out of ", total_graphs)
        R0, T, prob, size, AM, AF, JX,AM_t, AF_t, JX_t, AMAX, AFAX =sim.many_simulations_perc(100, G, t)
        
        T_list_f.append(t)
        synch_size.append(size)
        synch_prob.append(prob)
        synch_AM.append(AM)
        synch_JX.append(JX)
        synch_AF.append(AF)
        R0s_list.append(R0)
        AM_total.append(AM_t)
        AF_total.append(AF_t)
        JX_total.append(JX_t)
        AMAX_list.append(AMAX)
        AFAX_list.append(AFAX)
 
        graph_list.append(i)
        
    i = i +1


print("The average epidemic size across synch synthetic networks is: ", np.mean(synch_size))
print("THe average epidemic prob across synch synthetic networks is: ", np.mean(synch_prob)) 
print("THe average proportion of AMs infected across synch synthetic networks is: ", np.mean(synch_AM))    
print("THe average proportion of JXs infected across synch synthetic networks is: ", np.mean(synch_JX))  
print("THe average proportion of AFs infected across synch synthetic networks is: ", np.mean(synch_AF))  


epi_data = {"epi_size": synch_size,"Tr": T_list_f, "RO": R0s_list, "epi_prob": synch_prob, "graph": graph_list, "AM": synch_AM,  "JX": synch_JX, "AF": synch_AF, "AM_total": AM_total, "AF_total": AF_total, "JX_total": JX_total, "AM_AX": AMAX_list, "AF_AX": AFAX_list}
synch_epi_df = pd.DataFrame(epi_data)

synch_epi_df.to_csv("synch_epi_data_PC.csv")   
