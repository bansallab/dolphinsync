import random
import numpy as np
import math
import matplotlib.pyplot as plt
import networkx as nx
import csv
from numpy.random import uniform

def sort_dictionary(x):
# x is a dictionary
# x is sorted in descending order by value and a list of ordered keys is returned

    return sorted(x, key = x.get, reverse=True)


def probability_distribution_from_data(x):
# x is a data vector of positive values

	p = {}
	l = len(x)
	max_x = max(x)
	for i in range(0,max_x+1):
		p[i] = 0
	for xi in x:
		p[xi] = p[xi] + 1.0/l

	return p

def frequency_distribution_from_data(x):
# x is a data vector of positive values

	p = {}
	l = len(x)
	max_x = max(x)
	for i in range(0,max_x+1):
		p[i] = 0
	for xi in x:
		p[xi] = p[xi] + 1.0

	return p

def joint_probability_distribution_from_data(x,y):

	p = {}
	max_x = max(x)
	max_y = max(y)

	for i in range(0,max_x+1):
		for j in range(0,max_y+1):
			p[(i,j)] = 0

	total1 = 0
	total2 = 0
	for xi in x:
		for yj in y:
			p[(xi,yj)] = p[(xi,yj)] + 1
			total1 = total1 + 1

	for i in range(0,max_x+1):
		for j in range(0,max_y+1):
			p[(i,j)] = p[(i,j)]/(1.0*total1)

	return p

def joint_frequency_distribution_from_data(x,y):

	p = {}
	max_x = max(x)
	max_y = max(y)

	for i in range(0,max_x+1):
		for j in range(0,max_y+1):
			p[(i,j)] = 0

	total1 = 0
	total2 = 0
	for xi in x:
		for yj in y:
			p[(xi,yj)] = p[(xi,yj)] + 1

	return p

def entropy(px):

	e = 0

	for x in px.keys():
		if px[x] != 0:
			e = e + px[x]*np.log(1.0/px[x])

	return e

def normalized_entropy_of_degdist(px, N):

	e = entropy(px)
	return e/np.log(N-1)

def mutual_information(px, py, pxy):
	m = 0

	for y in py.keys():
		for x in px.keys():
			if py[y] != 0.0 and pxy[(x,y)] != 0.0:
				m = m + (pxy[(x,y)]*np.log(pxy[(x,y)]/py[y]) - px[x]*np.log(px[x]))

	return m

def variation_in_information(nx,ny,nxy, N):
	v = 0

	for x in nx.keys():
		for y in ny.keys():
			if nxy[(x,y)] > 0.0:
				v = v + ((nxy[(x,y)]/(1.0*N))*math.log((nx[x]*ny[y])/(nxy[(x,y)]**2)))

	v = v*(1.0/math.log(N))
	return v

def most_common(lst):
    return max(set(lst), key=lst.count)

def argmax_dict(dictt):
	key = []
	if dictt:
		m = max(dictt.values())
	return [key for key,val in dictt.iteritems() if val == m]

def remove_from_dict_key(dictt, remove_key):
	return {key: value for key, value in dictt.items() if key != remove_key}

def remove_from_dict_value(dictt, remove_value):
	return {key: value for key, value in dictt.items() if value != remove_value}

def read_dict_from_file(filestr):
	dictt = {}
	filetoread = open(filestr, 'r')
	csvreader = csv.reader(filetoread, delimiter=',')
	for row in csvreader:
		dictt[int(row[0])] = int(row[1])

	return dictt

def write_dict_to_file(dictt, filestr):
	filetowrite = open(filestr, 'w')
	csvwriter = csv.writer(filetowrite, delimiter=',')
	for key, val in dictt.items():
		csvwriter.writerow([key, val])

	return dictt

def draw_weighted_graph(G, node_colors, node_sizes,outstr):
# draws nodes in color list provided (node_colors) and edges based on weights

	pos=nx.random_layout(G) # positions for all nodes

	# draw nodes
	m = max(node_sizes.values())
	s = [(1000/m)*sz for sz in node_sizes.values()] # we want 1000 to be the largest node size
	nx.draw_networkx_nodes(G,pos,node_size=s, node_color=node_colors)



	# draw edges
	weights= [G[u][v]['weight'] for u,v in G.edges()]
	max_weight = max(weights)
	e = {}
	num_widths = 10
	for i in range(0,num_widths):
		e[i]=[(u,v) for (u,v,d) in G.edges(data=True) if (((1.0*d['weight'])/max_weight)>(i/(num_widths*1.0)) and ((1.0*d['weight'])/max_weight)<=((i+1.0)/num_widths))]
		nx.draw_networkx_edges(G,pos,edgelist=e[i], width=(i+1.0))

	plt.axis('off')
	plt.savefig(outstr) # save as png
	plt.show() # display

def slice_sampler(px, N = 1, x = None):
    """
    Provides samples from a user-defined distribution.

    slice_sampler(px, N = 1, x = None)

    Inputs:
    px = A discrete probability distribution.
    N  = Number of samples to return, default is 1
    x  = Optional list/array of observation values to return, where prob(x) = px.

    Outputs:
    If x=None (default) or if len(x) != len(px), it will return an array of integers
    between 0 and len(px)-1. If x is supplied, it will return the
    samples from x according to the distribution px.
    """
    values = np.zeros(N, dtype=int)
    samples = np.arange(len(px))
    px = np.array(px) / (1.*sum(px))
    u = uniform(0, max(px))
    for n in range(N):
        included = px>=u
        choice = random.sample(range(np.sum(included)), 1)[0]
        values[n] = samples[included][choice]
        u = uniform(0, px[included][choice])
    if x:
        if len(x) == len(px):
            x=np.array(x)
            values = x[values]
        else:
            print("px and x are different lengths. Returning index locations for px.")
    if N == 1:
        return values[0]
    return values

def dot_product(x,y):
	"""returns sum(x*y)"""
	return sum(list_product(x,y))

def dot_quotient(x,y):
	"""return sum(x/y)"""
	return sum(list_quotient(x,y))

def list_product(x,y):
	# returns (x*y)
	return ([i[0] * float(i[1]) for i in zip(x,y)])

def list_quotient(x,y):
	# return x/y
	return ([float(i[0])/float(i[1]) for i in zip(x,y)])

def flatten_matrix(matrix):
	""" Takes lists of lists and flattens to a combined list."""
	return [e for r in matrix.values() for e in r]

def deglist_to_degdist(deglist, start):
# returns degree distribution starting at degree = start
	n = len(deglist)
	return [float(list(deglist).count(i))/float(n) for i in range(start,max(deglist)+1)]

def prettyfloat(x):
# returns first three digits past decimal of a float
    return int(x*1000)/1000.0