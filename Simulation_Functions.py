#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 09:03:11 2021

@author: melissacollier
"""
import random as rnd             
import networkx as nx            
import statistics as st
import numpy as np
import os
import random


def get_networks_with_directory(directory):
    """
    This function returns a list of networks from a directory with multiple .graphml files

    """
    
    #first get files
    filelist = [filename for filename in sorted(os.listdir(os.path.abspath(directory))) if filename.endswith(".graphml")]
    networks = []
    os.chdir(directory)
    for file in filelist:
        network = nx.read_graphml(file)
        if len(network.nodes()) > 0:
            
            networks.append(network)
    return(networks)


def get_networks_with_directory_specifc(directory, contains):
    """
    This function returns a list of networks from a directory with multiple .graphml files
    Used if there are multiple network types in one directory and file names are 
    differentiated at the beginning
    """
    
    #first get files
    filelist = [filename for filename in sorted(os.listdir(os.path.abspath(directory))) if filename.endswith(".graphml")]
    networks = []
    os.chdir(directory)
    for file in filelist:
        if contains in file:
            network = nx.read_graphml(file)
            if len(network.nodes()) > 0:
            
                networks.append(network)
    return(networks)


######################################################################
def get_tau(G, R0, IP):
   
    """
    define the set parameters for the simulations:
   1. The average excess degree the number of edges connected to a neighbor, other than
      the edge the infection arrived along (essentially degree - 1). The average excess degree for the network is 
      the average of every node's excess degree.
   2. R0 is the number of new infections from one infectid
   3. T is the transmissibility of a pathogen in that network which is equal to R0/average excess degree
   4. gamma probability of recovery and is = 1/infectious period
   5. the return is tau, or the probability of infections and is derived from T = tau / (tau+gamma)
    """    
    
    nodes_names=list(G.nodes())
    #nodes= len(nodes_names)

    #######Find the average excess degree for the network to standardize R0#####
    
    
    x_deg_list = []
    for n in nodes_names:
        degree = G.degree(n)
        x = degree-1
        x_deg_list.append(x)

    
    ##For homogenous degree
    avg_x_deg = round(np.mean(x_deg_list),3)

    ### Now calculate tau
    T = R0/avg_x_deg 
    gamma = 1/IP #Recovery rate
    tau =(- T * gamma)/ (T - 1) #Transmission parameter
    
    return(tau, avg_x_deg)

def infected_neighbors(G, node, infected_list):
    """Calculates the number of infected neighbors for every susceptible node"""
    infected_deg = [x for x in list(G.neighbors(node)) if x in infected_list]   
    length = len(infected_deg)
    return length



def simulation_perc(Network,T):
# This function implements a percolation simulation to find R0
       
    net_size = Network.number_of_nodes()
    # Initialize variables for the list of infected and recovered individuals
    infected = []
    recovered = []

    ##################
    # Choose one node to infect (patient zero) so that outbreak can be seeded
    p_zero = rnd.choice(list(Network.nodes())) # Randomly choose one node from the network
    infected = [p_zero]                  # The node p_zero is now infected
    infected_count = 1
    R0_count = []  
  
    gen1_infected =[]
    while infected:
        
        infector = infected[0]      
        for neigh in list(Network.neighbors(infector)): # for all the nodes connected to (i.e. neighbors of) the infector 

            if neigh not in infected and neigh not in recovered: # check if this neighbor is susceptible

                # figure out if infector is successful at infecting neighbor "neigh"
                if rnd.random() < T:                     # if infector does infect neigh
                    gen1_infected.append(neigh) #append them to the gen1 infected list
                    infected_count = infected_count +1
        
        infected.remove(infector)
        recovered.append(infector)

######### For generation 2 ###############    
    infected = gen1_infected
    #print(infected)
    gen2_infected = []
    
    if len(infected) > 0:
        
        while infected:       
            num_infected = 0
            infector = infected[0]
            for neigh in list(Network.neighbors(infector)): # for all the nodes connected to (i.e. neighbors of) the infector 
                
                if neigh not in infected and neigh not in recovered and neigh not in gen2_infected: # check if this neighbor is susceptible

                # figure out if infector is successful at infecting neighbor "neigh"
                    if rnd.random() < T:                     # if infector does infect neigh
                        gen2_infected.append(neigh)
                        infected_count = infected_count +1
                        num_infected = num_infected +1
       
            R0_count.append(num_infected) # add the number they infected to the R0 list
            infected.remove(infector)
            recovered.append(infector)
        
        infected = gen2_infected
    else: infected = []
    
##### For gen 3 ################## 
    gen3_infected = []
    if len(infected) > 0:

        while infected:
            num_infected = 0
            infector = infected[0]
            for neigh in list(Network.neighbors(infector)): # for all the nodes connected to (i.e. neighbors of) the infector 
               
                if neigh not in infected and neigh not in recovered and neigh not in gen3_infected: # check if this neighbor is susceptible

                # figure out if infector is successful at infecting neighbor "neigh"
                    if rnd.random() < T:                     # if infector does infect neigh
                        gen3_infected.append(neigh)
                        infected_count = infected_count +1
                        num_infected = num_infected +1
            
            R0_count.append(num_infected) # add the number they infected to the R0 list
            infected.remove(infector)
            recovered.append(infector)

    
            
        infected = gen3_infected
    else: infected = []

    #FOR THE REMAINIG SIMULATION
    if len(infected) > 0:
        while infected:        
            infector = infected[0]
            for neigh in list(Network.neighbors(infector)): # for all the nodes connected to (i.e. neighbors of) the infector 
                
                if neigh not in infected and neigh not in recovered: # check if this neighbor is susceptible
#
                    if rnd.random() < T:                     
                        infected.append(neigh)
                        infected_count = infected_count +1        
#
            infected.remove(infector)
            recovered.append(infector)

    if len(R0_count) > 0:
        R0 = round(np.mean(R0_count),1)
    else: R0 = 0
    epi_size = infected_count/net_size
    # return R0, T value, and size of epidemic
    return R0, T, epi_size, recovered

def get_xdeg_with_deghet(G):
    
    #isolates = list(nx.isolates(G))
    #G.remove_nodes_from(isolates)
    #giant_cc = max(nx.connected_component_subgraphs(G), key=len)
    #G= giant_cc
    deg_list = [G.degree(node) for node in list(G.nodes())]
    avg_deg = round(np.mean(deg_list),3)
    var_deg = round(st.variance(deg_list),3)
    numer = var_deg + (avg_deg ** 2) - avg_deg
    xdeg = numer/avg_deg
    
    return xdeg

def simulation_perc_simple(Network,R0):
# This function implements a percolation simulation to find R0
       
    net_size = Network.number_of_nodes()
    # Initialize variables for the list of infected and recovered individuals
    
    xdeg = get_xdeg_with_deghet(Network)
    
    T = R0/xdeg
    ##################
    # Choose one node to infect (patient zero) so that outbreak can be seeded
    p_zero = rnd.choice(list(Network.nodes())) # Randomly choose one node from the network
    infected = [p_zero]                  # The node p_zero is now infected
    infected_count = 1 
    recovered = []
  
    while infected:
        
        infector = infected[0]      
        for neigh in list(Network.neighbors(infector)): # for all the nodes connected to (i.e. neighbors of) the infector 

            if neigh not in infected and neigh not in recovered: # check if this neighbor is susceptible

                # figure out if infector is successful at infecting neighbor "neigh"
                if rnd.random() < T:                     # if infector does infect neigh
                    infected.append(neigh) #append them to the gen1 infected list
                    infected_count = infected_count +1
        
        infected.remove(infector)
        recovered.append(infector)

    epi_size = infected_count/net_size
    # return R0, T value, and size of epidemic
    return T, epi_size, recovered



def many_simulations_perc(num_sims,G, T):
    """
    This function finds an average R0 value for a certain T value on a network
    It will only calculate R0 if the epidemic probability is 10%
    """   
    R0_list = []
    
    epidemic = []
    outbreak = []
    
    infected_AM = []
    infected_AF= []
    infected_JX = []
    
    AM_outof_total = []
    AF_outof_total = []
    JX_outof_total = []
    
    AM_AX_total = []
    AF_AX_total = []
    
   
    for x in range(num_sims):
        #print(x)
        R0, T, size, recovered = simulation_perc(G,T)
        
        if size >= 0.1: #only want epidemics so 10% or more of the network infected
            R0_list.append(R0)
            epidemic.append(size)
            
            
            AM = [] #rename or edit for each attribute type
            JX = []
            AF = []

            for node in G.nodes():
                if G.nodes[node]['Demo'] == 0: #change to relevent attribute type
                    AM.append(node)
                elif G.nodes[node]['Demo'] == 1:
                    AF.append(node)
                else: JX.append(node)
            AM_total = len(AM)
           # print("total AM: ", AM_total)
            JX_total = len(JX)
            #print("total JX: ", JX_total)
            AF_total = len(AF)
           # print("total AF", AF_total)
            
            ## Now get the total number infected for each attribute
            AM_inf = []
            JX_inf = []
            AF_inf = []

            
            for node in recovered:
                if G.nodes[node]['Demo'] == 0:
                    AM_inf.append(node)
                elif G.nodes[node]['Demo'] == 1:
                    AF_inf.append(node)
                else: JX_inf.append(node)
        
            AMinf_total = len(AM_inf)
            JXinf_total = len(JX_inf)
            AFinf_total = len(AF_inf)
            AXinf_total = AMinf_total +AFinf_total
            
            total_inf = AMinf_total + JXinf_total +AFinf_total
           
            AM_out_of_inf = AMinf_total /total_inf
            AF_out_of_inf = AFinf_total /total_inf
            JX_out_of_inf = JXinf_total /total_inf
            
            AM_out_of_infAX = AMinf_total/AXinf_total
            AF_out_of_infAX = AFinf_total/AXinf_total
            
            AM_outof_total.append(AM_out_of_inf)
            AF_outof_total.append(AF_out_of_inf)
            JX_outof_total.append(JX_out_of_inf)
            
            AM_AX_total.append(AM_out_of_infAX)
            AF_AX_total.append(AF_out_of_infAX)
            
            
            propAM = round(AMinf_total/AM_total, 2)
            propJX = round(JXinf_total/JX_total, 2)
            propAF = round(AFinf_total/AF_total, 2)
           
            infected_AM.append(propAM)
            infected_JX.append(propJX)
            infected_AF.append(propAF)
        ###### END ATTRIBUTE CODE ##################
        #############################################################3
        else: 
            outbreak.append(size)
            
                 
    if len(epidemic) == 0:
    
        avg_epi_size = 0
        prob_epidemic = 0
        avg_prop_AM = 0
        avg_prop_AF = 0
        avg_prop_JX =0
        R0 = 0
        AMoototal = 0
        AFoototal = 0
        JXoototal = 0
        AM_AX = 0
        AF_AX = 0
     
        
    else: 
        avg_epi_size = round(np.mean(epidemic),2)
        
        num_epidemics = len(epidemic)
        prob_epidemic = round((num_epidemics/num_sims),2)
        
        ########## Comment out if no attributes ###################
        avg_prop_AM = round(np.mean(infected_AM),2)
        avg_prop_JX = round(np.mean(infected_JX),2)
        avg_prop_AF = round(np.mean(infected_AF),2)
        
        AMoototal = round(np.mean(AM_outof_total), 2)
        AFoototal = round(np.mean(AF_outof_total), 2)
        JXoototal = round(np.mean(JX_outof_total), 2)
        
        AM_AX = round(np.mean(AM_AX_total), 2)
        AF_AX = round(np.mean(AF_AX_total), 2)

        
        ###########################################################

        R0 = round(np.mean(R0_list),1)  
        

    return R0, T, prob_epidemic, avg_epi_size, avg_prop_AM, avg_prop_AF,avg_prop_JX, AMoototal, AFoototal, JXoototal, AM_AX, AF_AX
    



def chain_binomial_simulation(G,tau, gamma, trange):
    """
    This function is a typical SIR simulation using the chain binomial method.
    """
# Initialize variables for the list of infected and recovered individuals
    infected = []
    recovered = []
    net_size = len(G.nodes())
    ##################
    # Choose one node to infect (patient zero)
    
    p_zero = list(np.random.choice(list(G.nodes()), 1)) # Randomly choose one node from the network
    
    infected = list(p_zero)                  # The node p_zero is now infected
    infected_nodes = list(p_zero)
    ##################
    # Run the simulation to simulate disease spread in network
        
    for t in range(trange): ##for every day of the epidemic
        
        for node in list((G.nodes)):

            if node not in infected and node not in recovered:
                
                prob = 1- np.exp(-tau[0]*infected_neighbors(G, node, infected))
                
                if rnd.random() < prob:
                    infected.append(node)
                    infected_nodes.append(node)
        
        for infector in infected:
            if rnd.random() < gamma:
                infected.remove(infector)
                recovered.append(infector)

        ##################
        #num_infected = len(infected)/net_size
    size = len(recovered)/net_size
 
    
    return size, infected_nodes

def many_simulations(num_sims, G, tau, IP, days):
    """
    This function simulates many epidemics on one network. It can also show
    a proportion of infected individs by attribute (comment out section if 
                                                     no node attributes)
    """   
    epidemic = []
    outbreak = []
    
    infected_AM = []
    infected_JX= []
    infected_AF = []
    
    gamma = 1/IP
    
    
   
    for x in range(num_sims):
        print(x)
        size, infected_list = chain_binomial_simulation(G, tau, gamma, days)
        if size > 0.1:
            epidemic.append(size)
        
        
        ######Comment out section below if no attributes################
        #################################################################
        ####First find total of each attribute type in network#######
            AM = [] #rename or edit for each attribute type
            JX = []
            AF = []
            
            for node in G.nodes():
                if G.nodes[node]['Demo'] == 0: #change to relevent attribute type
                    AM.append(node)
                if G.nodes[node]['Demo'] == 1:
                    AF.append(node)
                else: JX.append(node)
            AM_total = len(AM)
            JX_total = len(JX)
            AF_total = len(AF)
        
        ## Now get the total number infected for each attribute
            AM_inf = []
            JX_inf = []
            AF_inf = []
            
        
            for node in infected_list:
                if G.nodes[node]['Demo'] == 0:
                    AM_inf.append(node)
                if G.nodes[node]['Demo'] == 1:
                    AF_inf.append(node)
                else: JX_inf.append(node)
        
            AMinf_total = len(AM_inf)
            JXinf_total = len(JX_inf)
            AFinf_total = len(AF_inf)
            
            propAM = round(AMinf_total/AM_total, 2)
            propJX = round(JXinf_total/JX_total, 2)
            propAF = round(AFinf_total/AF_total, 2)
            
            infected_AM.append(propAM)
            infected_JX.append(propJX)
            infected_AF.append(propAF)
            
        ###### END ATTRIBUTE CODE ##################
        #############################################################3
        else: outbreak.append(size)
             
    if len(epidemic) == 0:
        avg_epi_size = 0
    else: avg_epi_size = round(np.mean(epidemic),1)
    
    num_epidemics = len(epidemic)
    prob_epidemic = round((num_epidemics/num_sims),2)
    
    ########## Comment out if no attributes ###################
    avg_prop_AM = round(np.mean(infected_AM),2)
    avg_prop_JX = round(np.mean(infected_JX),2)
    avg_prop_AF = round(np.mean(infected_AF),2)

    ###########################################################
    ##### remove attribute from return if not using
    
    
    
    return avg_epi_size, prob_epidemic, avg_prop_AM, avg_prop_JX, avg_prop_AF 

    
      
    