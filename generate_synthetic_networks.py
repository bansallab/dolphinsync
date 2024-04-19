#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 13:46:36 2023

@author: melissacollier

This code can generate networks from a specified degree distriution, mixing matrix and demographic distribution
As downloaded this code will generate 25 graphs using PC data into a folder called "generated_graphs". 
"""

import networkx as nx
import random
import numpy as np
import time
import general_tools
import pretty_print
from scipy import special
import pandas as pd
from functools import reduce
 
# Function to get unique values
 
 
def unique(list1):
 
    # Print directly by using * symbol
    ans = reduce(lambda re, x: re+[x] if x not in re else re, list1, [])
    print(ans)


def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)

def estimate_num_edges_M(N, node_type_distribution, degree_dist_type):
    num_types = len(degree_dist_type)
    node_type_frequency = [round(N*p) for p in node_type_distribution]
    average_degree_by_type = [general_tools.dot_product(range(1,1+len(degree_dist_type[i])),degree_dist_type[i]) for i in range(0,num_types)]
    M = general_tools.dot_product(node_type_frequency, average_degree_by_type)/2.0

    return round(M)

def sample_edges(M, mixing_matrix):

	# sample M edge types from mixing matrix
	flat_mixing = general_tools.flatten_matrix(mixing_matrix)
	edge_types = general_tools.slice_sampler(flat_mixing, N=int(M))

	# take edges types and find origin destination node types
	index_lookup = {}
	index_lookup = [(i,j) for i in range(0, len(mixing_matrix)) for j in range(0, len(mixing_matrix))]
	sampled_edges = [index_lookup[e] for e in edge_types]

	return sampled_edges


def count_edges_of_type_i(list_of_edges):
	node_types_in_edges = [i for t in list_of_edges for i in list(t)]
	#node_types_in_edges = [t1 for t1,t2 in list_of_edges] # count start of edge only
	m_i = general_tools.frequency_distribution_from_data(node_types_in_edges)
	return m_i.values()


def create_nodes(G, num_types, num_nodes_by_type):
    """ Add to the graph G, nodes of each type.  """
    node_id = 0
    print(num_nodes_by_type)
    print(num_types)
    for node_type in range(0, num_types):
        for t in range(0, num_nodes_by_type[node_type]):
            G.add_node(node_id, node_type=node_type)
            node_id += 1
    return G

def assign_degrees(G, degree_dist_type, num_types, numedges_type):
	""" Assign degrees by type to G. Ensure number of edges is adhered to.  """

	# assign node degree
	attributes = nx.get_node_attributes(G, 'node_type')
	for n_type in range(0, num_types):
        
		nodes_of_type = [k for k, a in attributes.items() if a==n_type]
		max_deg = len(degree_dist_type[n_type])
		degree_list = general_tools.slice_sampler(degree_dist_type[n_type], N = len(nodes_of_type), x = range(1, max_deg+1))
		node_degrees = {x[0]:x[1] for x in zip(nodes_of_type, degree_list)}
		nx.set_node_attributes(G, name = 'degree', values = node_degrees)

		sum_of_degrees_type = sum([G.nodes[k]['degree'] for k in nodes_of_type])
		while sum_of_degrees_type != numedges_type[n_type]: # sum_of_degrees_expected = numedges_type[n_type]
			node = random.choice(nodes_of_type)
			deg= general_tools.slice_sampler(degree_dist_type[n_type], N = 1, x = range(1, max_deg+1))

			# need to make sure newly assigned degrees are helping move the sum closer to the desired sum
			cond1 = ((sum_of_degrees_type - numedges_type[n_type]) > 0)
			cond2 = ((G.nodes[node]['degree'] - deg) > 0)
			while int(cond1)+int(cond2) == 1:
				node = random.choice(nodes_of_type)
				deg= general_tools.slice_sampler(degree_dist_type[n_type], N = 1, x = range(1, max_deg+1))
				cond2 = ((G.nodes[node]['degree'] - deg) >= 0)
			G.nodes[node]['degree'] = deg
			sum_of_degrees_type = sum([G.nodes[k]['degree'] for k in nodes_of_type])

	return G

def find_edge_list(candidate_edge_list, G):

	edge_list = []
	for n1,n2 in candidate_edge_list:
		if n1 != n2 and n2 not in G.neighbors(n1):
			edge_list.append((n1,n2))

	return edge_list

def assign_edges(G,mxmatrix, num_types):
# code originally written to sort nodes by degree (ie use Havel-Hakimi). Not doing that anymore

	# make copies of nodes of each type so can connect them
	node_copies = {} # contains nodes for each type sorted by degree
	for ntype in range(0,num_types):
		nodes_type_degree = [(n,G.nodes[n]['degree']) for n in G.nodes() if G.nodes[n]['node_type']==ntype]
		#nodes_degree_sorted = sorted(nodes_type_degree,key=lambda x: x[1])
		#copy_lists = [deg*[node] for node, deg in nodes_degree_sorted]
		copy_lists = [deg*[node] for node, deg in nodes_type_degree]
		node_copies[ntype] = [i for xx in copy_lists for i in xx] # collapse from list of lists to list

	# link nodes by type
	for ntype1 in range(0,num_types):
		for ntype2 in range(0,num_types):
			while mxmatrix[ntype1][ntype2] > 0:
				node2 = random.choice(node_copies[ntype2]) # get a random item for node2
				node_copies[ntype2].remove(node2) # and remove it from the list
				node1 = random.choice(node_copies[ntype1]) #  get a random item for node1
				node_copies[ntype1].remove(node1) # and remove it from the list
				while node1 == node2 or node2 in G.neighbors(node1):
				# NOTE: THIS LOOP GETS STUCK WHEN THERE ARE ONLY A COUPLE OF ITEMS LEFT IN THE LIST
					node_copies[ntype2].insert(0,node2) # put the removed item back
					node_copies[ntype1].insert(0,node1) # put the removed item back
					node2 = random.choice(node_copies[ntype2]) # and get a random item
					node_copies[ntype2].remove(node2) # and remove it from the list
					node1 = random.choice(node_copies[ntype1]) # and get a random item
					node_copies[ntype1].remove(node1) # and remove it from the list
				mxmatrix[ntype1][ntype2] = mxmatrix[ntype1][ntype2]-1
				G.add_edge(node1, node2)

	return G

def connected_double_edge_swap_of_type(G, type1, type2):
	candidate_edges = [(i,j) for i,j in G.edges() if check_type(G, i, j, type1, type2)]
	edge1 = random.choice(candidate_edges)
	edge2 = random.choice(candidate_edges)

	G.remove_edges_from([edge1, edge2])
	type1_edge1 = G.node[edge1[0]]['node_type']
	type1_edge2 = G.node[edge2[0]]['node_type']
	if type1_edge1 != type1_edge2:
		G.add_edges_from([(edge1[0],edge2[0]),(edge1[1],edge2[1])])
	else:
		G.add_edges_from([(edge1[0],edge2[1]),(edge1[1],edge2[0])])

	return G

def randomize_graph(G, num_types):
# randomize a network using double-edged swaps
# note this swapping algorithm does not check for connectivity

	size = G.size() # number of edges in graph
	its = 5*size
	for type1 in range(0,num_types):
		for type2 in range(0, num_types):
			for counter in range(0,its):
				G = connected_double_edge_swap_of_type(G, type1, type2)

	return G

def check_type(c, i, j, type1, type2):
	return (c.nodes[i]['node_type'], c.nodes[j]['node_type']) == (type1,type2) or\
	(c.nodes[i]['node_type'], c.nodes[j]['node_type']) == (type2,type1)


def get_candidate_edges_of_type(c, type1, type2):
	""" Given a component, return an edge where both nodes have degree 1 and
	types match. """

	at_least_degree_two = [(i,j) for i, j in c.edges() if c.degree(i) > 1 and c.degree(j) > 1]
	candidate_edges = [(i,j) for i,j in at_least_degree_two if check_type(c, i, j, type1, type2)]
	if candidate_edges:
	   return random.choice(candidate_edges)
	else:
	   return None

def get_candidate_edge(component):
	""" Given a component, return an edge where either node has degree 1 (so
	component won't get disconnected when we re-wire the graph). If there's only
	one edge in the component, return that. """

	edges = list(component.edges())
	if len(edges) > 1:
		at_least_degree_two = [(i,j) for i,j in edges if component.degree(i) > 1 or component.degree(j) > 1]
		edge = random.choice(at_least_degree_two)
		return edge
	else:
		return(edges[0])


def connect_graph(G):
	""" Check if G is disconnected and connect if necessary using modified
	Taylor's algorithm """

	cc = list(connected_component_subgraphs(G))
	component_count = len(cc)
	#sz = [c.number_of_nodes() for c in cc]
	#print sz

	while component_count > 1:   #while G is not connected, reduce number of components

    #print(component_count)

		# pick a random edge in the second largest component (cc[1])
		edge2 = get_candidate_edge(cc[1])
		type1 = G.nodes[edge2[0]]['node_type']
		type2 = G.nodes[edge2[1]]['node_type']
		print(type1, type2)

		# pick a random edge in the largest component (cc[0]) where both nodes
		# have degree > 1 and the types match those of edge2
		edge1 = get_candidate_edges_of_type(cc[0], type1, type2)
		while edge1 is None:
		  edge2 = get_candidate_edge(cc[1])
		  type1 = G.nodes[edge2[0]]['node_type']
		  type2 = G.nodes[edge2[1]]['node_type']
		  edge1 = get_candidate_edges_of_type(cc[0], type1, type2)
		  #print edge1, edge2

		# swap connections between edge1 and edge2
		#  to attempt to connect the two components
		G.remove_edges_from([edge1, edge2])
		type1_edge1 = G.nodes[edge1[0]]['node_type']
		type1_edge2 = G.nodes[edge2[0]]['node_type']
		if type1_edge1 != type1_edge2:
			G.add_edges_from([(edge1[0],edge2[0]),(edge1[1],edge2[1])])
		else:
			G.add_edges_from([(edge1[0],edge2[1]),(edge1[1],edge2[0])])

		cc = list(connected_component_subgraphs(G))
		component_count = len(cc)
	return G


def cleanup_graph(G):
# remove self loops (there shouldn't be any multiedges in Graph())
# CODE NEEDS TO BE IMPROVED

	# replace self loops with like-type edges
	cc = 0
	while len(list(nx.selfloop_edges(G))) > 1 and cc < 10*len(list(nx.selfloop_edges(G))):
		print("selfloops:", G.selfloop_edges())
		edge1 = random.choice(G.selfloop_edges())
		edge2 = random.choice(G.selfloop_edges())
		c = 0
		while (edge1[0] == edge2[0] or G.node[edge1[0]]['node_type'] != G.node[edge2[0]]['node_type']) and c < 100:
			edge2 = random.choice(G.selfloop_edges())
			c = c + 1
		G.remove_edge(edge1[0],edge1[1])
		if edge2 in G:
			G.remove_edge(edge2[0],edge2[1])
		G.add_edge(edge1[0], edge2[0])
		cc = cc + 1

	# if rewiring wasn't successful just remove the selfloops
	for e1,e2 in list(nx.selfloop_edges(G)):
		G.remove_edge(e1,e2)

	return G

def check_graph_characteristics(G, degree_dist_type, mixing_matrix, num_types):

	degdist = {}
	mxmatrix = {}
	nodetypes = nx.get_node_attributes(G, 'node_type')
	avg_deg_desired = []
	avg_deg_actual = []
	for i in range(0, num_types):
		nodes_of_type = [node for node, ntype in nodetypes.items() if ntype==i]
		degdist[i] = general_tools.deglist_to_degdist(dict(G.degree(nodes_of_type)).values(), 1)
		avg_deg_desired.append(general_tools.dot_product(range(1,len(degree_dist_type[i])),degree_dist_type[i]))
		avg_deg_actual.append(general_tools.dot_product(range(1,len(degdist[i])),degdist[i]))

	print("expected degree dist: \n", pretty_print.print_dictionary(degree_dist_type))
	print("actual degree dist: \n", pretty_print.print_dictionary(degdist))

	print("expected avgdegree: \n", pretty_print.print_list(avg_deg_desired))
	print("actual avg degree: \n", pretty_print.print_list(avg_deg_actual))
	print("total network avg degree: \n", sum(dict(G.degree()).values())/float(G.number_of_nodes()))

	num_edges = G.number_of_edges()
	for n1,n2 in G.edges():
		type1 = G.nodes[n1]['node_type']
		type2 = G.nodes[n2]['node_type']
		if type1 not in mxmatrix:
			mxmatrix[type1] = {}
		if type2 not in mxmatrix[type1]:
			mxmatrix[type1][type2] = 0
		mxmatrix[type1][type2] = mxmatrix[type1][type2] + 1.0/num_edges

	print("expected mixing matrix: \n",pretty_print.print_dictionary(mixing_matrix))
	print("acutal mixing matrix: \n", pretty_print.print_dictionary(mxmatrix))


def generate_graph(N, node_type_distribution, degree_dist_type, mixing_matrix):#, weight_distribution_typetype):
	# N: number of nodes in network
	# node_type_distribution: specifies distribution of node types (e.g. calf, juvenile, adult, etc.)
	# degree_dist_type: each row specifies degree distribution for each node type, starts at degree = 1
	# mixing matrix: each row specifies proportion of edges to each node type, whole matrix sums up to 1
	# weight_distribution_typetype: provides distribution of edge weights for each pair of

	start_time = time.time()

	num_types = len(degree_dist_type)

	# initialize graph
	G = nx.Graph()

	# figure out how many edges are needed
	M = estimate_num_edges_M(N, node_type_distribution, degree_dist_type)
	print("number of edges: ", M,"\n")
	print("mixing_matrix:\n")
	pretty_print.print_dictionary(mixing_matrix)

	# sample M edges from mixing_matrix (returns a list of node-type tuples)
	sampled_edges = sample_edges(M, mixing_matrix)

	# calculated mxmatrix
	mxmatrix = {k:[0]*num_types for k in range(0,num_types)}
	for type1,type2 in sampled_edges:
		mxmatrix[type1][type2] = mxmatrix[type1][type2] + 1.0
    
	print("mxing_matrix calc from sampled edges:\n", pretty_print.print_dictionary(mxmatrix),"\n", type(mxmatrix))
    

	# for each nodetype, figure out m_i= number of edges that end in type i
	m_i = list(count_edges_of_type_i(sampled_edges))
    
	print("number of edges of each type: ", m_i,"\n")
	print("number of total_edges: ", sum(m_i),"\n")

	# find zi = avgdegree by type and ni = mi/zi
	average_degree_by_type = [general_tools.dot_product(range(1,1+len(degree_dist_type[i])),degree_dist_type[i]) for i in range(0,num_types)]
	n_i = general_tools.list_quotient(m_i, average_degree_by_type)
	n_i = [n for n in n_i] # need to correct for the fact that there are m_i*2 edges of each type
	n_i = list(map(int, map(round,n_i)))
	print("average degree by type: \n", pretty_print.print_list(average_degree_by_type))
	print("number of total nodes: \n", sum(n_i),"\n")
	print("expected node type distribution: \n", pretty_print.print_list(node_type_distribution))
	print("actual node type distribution: \n", pretty_print.print_list(general_tools.list_quotient(n_i, [sum(n_i)]*len(n_i))))

	# create nodes by type based on n_i
	print("creating nodes","\n")
	G = create_nodes(G, num_types, n_i)

	# for each nodetype, draw degrees from appropriate degree_dist
	print("assigning degrees","\n")
	G = assign_degrees(G, degree_dist_type, num_types, m_i)

	# assign edges a node based on degree
	print("assigning edges","\n")
	G = assign_edges(G, mxmatrix, num_types)

	# randomizing graph edges (by type). Only need if using Havel-hakimi type connections (which not doing)
	#print "randomizing edges"
	#G = randomize_graph(G, num_types)

	# make sure graph is connected
	print("connecting, cleaningup graph","\n")
	#print "isolated nodes: ", len([i for i,d in G.degree().items() if d == 0])
	G = connect_graph(G)
	G = cleanup_graph(G)

	print("Checking graph statistics","\n")
	check_graph_characteristics(G, degree_dist_type, mixing_matrix, num_types)

	print("Time spent: ", (time.time()-start_time)/60.0, " minutes","\n")

	# add weights on all graph edges

	return G



def generate_mixing_matrix(avg_list, sd_list):
  # takes an input of the average and std of a mixing matrix e.g. avg_list = [[0.25], [0.08,0.17], [0.12,0.21,0.17]]
  ##                                                              sd_list= [[0.12], [0.03, 0.025], [0.1, 0.04, 0.03]]
  
    num_rows = len(avg_list) 
    print(num_rows)
    ## Give you the number of rows for the mixing matrix

    matrix_rows = []

    for x in range(0, num_rows):
       
      
       avg = avg_list[x]
       sd = sd_list[x]
       i = 0
       row = []
       
       ## This for loop will generate values for each row in the mixing matrix
       ### And make sure that the values are = 1
       for mu in avg:
           
           std = sd[i]
           
           
           min_value = mu - std
           max_value = mu+std
           
           min_in_standard_domain = (min_value - mu) / std  
           max_in_standard_domain = (max_value - mu) / std
           
           if min_in_standard_domain <= 0:
               min_in_standard_domain == 0
           
           min_in_erf_domain = special.erf(min_in_standard_domain)
           max_in_erf_domain = special.erf(max_in_standard_domain)

           random_uniform_data = np.random.uniform(min_in_erf_domain, max_in_erf_domain, 1)
           random_gaussianized_data = (special.erfinv(random_uniform_data) * std) + mu
           x = abs(round(np.mean(random_gaussianized_data), 2))
           row.append(x)
               
       matrix_rows.append(row)
    
    
    

   
    return(matrix_rows)

             
    

def generate_degree_distribution(df, num_demos): 
    
##### Uses an exported fitted dsitribution from R
#### Makes sure the columnsa are in the same order you provide your mixing matrix data in
##### This order is AX = 0, AFNN = 1, AFNNN = 2, JX = 3
    
    dist_list = []
    
    max_list = []
    for demo in range(num_demos):
        m = df[df.columns[demo]].max()
        max_list.append(m)
    
    largest_deg = max(max_list)
    
    for demo in range(num_demos):
    
       d = df[df.columns[demo]].values.tolist()
       unique_degrees = list(set(d))
       print(unique_degrees)
       
       deg_list_demo = []
       
       for deg in range(largest_deg+3):
         
           
           if deg not in unique_degrees:
               deg_list_demo.append(0)
           else:
               occ = d.count(int(deg))/1000
              
               deg_list_demo.append(occ)

        
       dist_list.append(deg_list_demo)
    
    return(dist_list)

def specify_data_new(num_types, degree_data, mixing_data):
# updated Oct 2016

    node_type_dist = []
    degree_dist = {}
    mixing_matrix = {k:[0]*num_types for k in range(0,num_types)}

# specify size of network (number of nodes)
    N = 2000
# specify proportion of population belonging to each demographic groups as ordered below
    node_type_dist = [0.42, 0.44, 0.14] # AM, AF, JX

# For Degree Distribution
   # For Degree Distribution
    degree_dist = {}  
    degree_dist[0]= degree_data.iloc[:, 1].tolist()
    degree_dist[1]= degree_data.iloc[:, 2].tolist()
    degree_dist[2]= degree_data.iloc[:, 3].tolist()
# For mixing matrix
    
    mixing_values = mixing_data['PC_Mean_Ratio'].tolist()
    sd_values = mixing_data['PC_SD_Ratio'].tolist()
    mixing_matrix = {}# initialize the variable as a dictionary (where each value will be a list)
    avg_list = [[mixing_values[0], mixing_values[1], mixing_values[2]], [mixing_values[3], mixing_values[4], mixing_values[5]], [mixing_values[6], mixing_values[7], mixing_values[8]]]
    sd_list= [[sd_values[0],sd_values[1], sd_values[2]], [sd_values[3], sd_values[4], sd_values[5]], [sd_values[6], sd_values[7], sd_values[8]]]
    mixing_matrix[0], mixing_matrix[1],mixing_matrix[2]= generate_mixing_matrix(avg_list, sd_list)

    return N, node_type_dist, degree_dist, mixing_matrix

def output_graph_and_statistics(G, fileext, graphnum):

	# output edgelist with delimiter comma and no edge data
	nx.write_edgelist(G, 'generated_graphs/edgelist_'+fileext+'_'+str(graphnum)+'.txt', data=False)

	# ouput nodetype dictionary
	pretty_print.print_dictionary_to_file(nx.get_node_attributes(G,'node_type'), "generated_graphs/nodetypes_"+fileext+"_"+str(graphnum)+".txt")

	# output degree dist by type
	degdist = {}
	nodetypes = nx.get_node_attributes(G, 'node_type')
	num_types = len(set(nodetypes.values()))
	for i in range(0, num_types):
		nodes_of_type = [node for node, ntype in nodetypes.items() if ntype==i]
		degdist[i] = general_tools.deglist_to_degdist(dict(G.degree(nodes_of_type)).values(), 1)
	pretty_print.print_dictionary_to_file(degdist, "generated_graphs/type_degreedist_"+fileext+"_"+str(graphnum)+".txt")

	# output degre dist -- total
	p = general_tools.probability_distribution_from_data(dict(G.degree()).values())
	pretty_print.print_dictionary_to_file(p, "generated_graphs/total_degreedist_"+fileext+"_"+str(graphnum)+".txt")

	# ouptut mixing matrix
	# --> compute mixing matrix and then print
	mxmatrix = {}
	num_edges = G.number_of_edges()
	for n1,n2 in G.edges():
		type1 = G.nodes[n1]['node_type']
		type2 = G.nodes[n2]['node_type']
		if type1 not in mxmatrix:
			mxmatrix[type1] = {}
		if type2 not in mxmatrix[type1]:
			mxmatrix[type1][type2] = 0
		mxmatrix[type1][type2] = mxmatrix[type1][type2] + 1.0/num_edges
	pretty_print.print_dictionary_to_file(mxmatrix, "generated_graphs/mixing_matrix_"+fileext+"_"+str(graphnum)+".txt")


###############################################################
if __name__ == "__main__":

    num_graphs = 25 # specificy number of synthetic networks to create
    num_demog_groups = 3 # specify number of demographic groups you want network to have
	
    for graphnum in range(1,num_graphs+1):

        print("Generating graph number: ", graphnum)
        deg_data = pd.read_csv('PC_degree_distributions.csv')
        mm_data = pd.read_csv('Mixing_Matrix_data_PC.csv')

    # you'll have to go into this function and 
        [N , node_type_dist, degree_dist_type, mixing_matrix] = specify_data_new(num_demog_groups, deg_data,mm_data)
        G_before96 = generate_graph(N, node_type_dist, degree_dist_type, mixing_matrix)
        output_graph_and_statistics(G_before96, 'PC_network',graphnum)
        #output_graph_and_statistics(G_before96, 'Synchrony',graphnum)

