#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 14:33:23 2023

@author: melissacollier
"""

"""
To add attributes from a csv file to a graphml the following algorithim will be used:
    1. Import network
    2. Import attribute file
    3. Replace NA values with a physical blank values
    4. Create a dictionary for every column/node combo
    5. Add each dictionary as a node attribute using networkx
    6. Export the new graphml file with attributes.
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
import networkx as nx

dirct = str(Path.cwd())    
directory = sorted(os.listdir(os.path.abspath(dirct +"/generated_graphs")))

edgelists= []
node_attr = []
    
for file in directory:
    print(file)
    if file.startswith('edgelist'):
        edgelists.append(file)
       
            
    if file.startswith('nodetypes'):
        node_attr.append(file)


#print(edgelists)
#print(node_attr)


i = 0
for file in edgelists:
    os.chdir(dirct+"/generated_graphs")
    filenet = str(i)
    #net = str(net)
    G = nx.read_edgelist(file)
    attributes = pd.read_csv(node_attr[i], sep = ",", header =None )
    attributes.columns = ["Node", "Demo"]
    attributes["Node"] = attributes.Node.astype(str)
    att_dict = attributes.set_index('Node').T.to_dict()
    nx.set_node_attributes(G, att_dict)
    nx.write_graphml(G, "PC_Network"+ filenet + ".graphml" )
    i = i+1
   

   
   
   
   