# dolphinsync

# This repository contains data and code for "Breathing in sync: how a social behavior structures respiratory epidemic risk in bottlenose dolphins"

## Data files:</br>

1. **Clean_Degree_Data_PCDP.csv** contains focal follow data from PC. Column descriptions: </br>
   _Follow_ID_ = unique ID for the focal follow conducted</br>
   _Dolphin-ID_ = unique ID of the focal dolphin</br>
   _Demo_ = Demographic group: AM = Adult Male, AF= Adult Female, AFNN= Mother and neonate calf pair, AFNNN = mother and nonneonate calf pair, JX = juvenile</br>
   _Demo_Confidenc_e = Confidence in demo assignment based on Appdendix 1 in paper</br>
   _Degree_ = total sync degree during follow</br>
   _AM_deg_ = total degree during follow that was an AM</br>
   _AF_deg_ = total degree during follow that was an AF</br>
   _JX_deg_ = total degree during follow that was a JX</br>
   _low_AM_deg_ = total AM_deg that was low confidence demo assignment</br>
   _low_AF_deg_ = total AF_deg that was low confidence demo assignment</br>
   _low_JX_deg_ = total JX_deg that was low confidence demo assignment</br>
   _Follow_length_ = length of follow in minutes</br>
   _Date_ = date of follow</br>
   _Number_Syncs_ = number of unique synchronized breathing events that occurred in the follow</br>
   _Group_Size_ = total individuals present during the focal follow</br>
   _Year_ = year of follow</br>
2. **Clean_Degree_Data_SB.csv** contains focal follow data from SB. Column descriptions:</br>
   _Follow_ID_ = unique ID for the focal follow conducted</br>
   _Dolphin-ID_ = unique ID of the focal dolphin</br>
   _Demo_ = Demographic group: AM = Adult Male, AF= Adult Female, AFNN= Mother and neonate calf pair, AFNNN = mother and nonneonate calf pair, JX = juvenile</br>
   _Degree_ = total sync degree during follow</br>
   _AM_deg_ = total degree during follow that was an AM</br>
   _AF_deg_ = total degree during follow that was an AF</br>
   _JX_deg_ = total degree during follow that was a JX</br>
   _Follow_length_ = length of follow in minutes</br>
   _Date_ = date of the follow </br>
   _Number_Syncs_ = number of unique synchronized breathing events that occurred in the follow </br>
   _Group_Size_ = total individuals present during the focal follow</br>
   _Year_ = year of the follow </br>
3. **Stranding_Data_background_removed_by_Adult_Sex.csv** contains the stranding records during the UME with the background strandings removed for each adult sex class and state based on our methods</br>
4. **Stranding_Data_background_removed_by_Age.csv** contains the stranding records during the UME with the background strandings removed for each adult sex class and state based on our methods.</br>
5. **PC_degree_distributions.csv** and **SB_degree_distributions.csv** are generated using code file 1. Each row represents a degree (e.g. row 1 = degree 0, row 2 = degree of 1, etc) and each column (AM, AF, JX) shows the proportion of that demographic that would have that degree over an average DMV infectious period.</br>
6. **Mixing_Matrix_data_PC.csv** and  **Mixing_Matrix_data_SB.csv** are generated using code file 2. Column descriptions:</br>
   _Comb_: demo combination based on focal:contact</br>
   _Mean_: average degree for that combination across focal follows</br>
   _Lower_: the lower degree estimate of that combination</br>
   _Upper_: the upper estimate of that combination</br>
   _Mean_Ratio_: The proportion of all observed syncs that are of that combination acorss all follows </br>
   _SD-Ratio_: the std of that proportion </br>
   _empirical_lower_ratio/upper_ratio_: the upper and lower bounds of that proportion </br>
7. **generated_graphs**: a sub directory with all network files for 25 networks created based on PC data. All files in this sub directory were generated using code file 3 </br></br>

## Code Files:</br>

1. **Estimate_degree_distribution_for_infectious_period.R**: R code to generate a degree distribution for AM, AF and JX dolphins using data from PC and SB (data files 1 and 2) </br>
2. **Generate_Mixing_Matrices.R**: Code to generate a mixing matrix across demographic groups using the PC and SB data files.</br>
3. **generate_synthetic_networks.py**: Code to generate synthetic networks based on emprical degree distributions and mixing matrices generated code files 1 and 2 (i.e. data files 5 and 6). Requires the user make a subflder in their working directory called "generated_networks". Also requires code files **general_tools.py** and **pretty_print.py**</br>
4. **Create_graphml_networks_for_disease_Simulations.py**: Code to create a graphml file for all generated networks to be used for disease simulations.</br>
5. **Disease_Simulation_Code.py**: Code to run disease simulations on each generated network and write all model results into a csv file. Requires the code file: **Simulation_Functions.py**</br>


