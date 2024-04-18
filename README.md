# dolphinsync

This repository contains data and code for "Breathing in sync: how a social behavior structures respiratory epidemic risk in bottlenose dolphins"

Data files:

1. Clean_Degree_Data_PCDP.csv contains focal follow data from PC. Column descriptions:
   Follow_ID = unique ID for the focal follow conducted
   Dolphin-ID = unique ID of the focal dolphin
   Demo = Demographic group: AM = Adult Male, AF= Adult Female, AFNN= Mother and neonate calf pair, AFNNN = mother and nonneonate calf pair, JX = juvenile
   Demo_Confidence = Confidence in demo assignment based on Appdendix 1 in paper
   Degree = total sync degree during follow
   AM_deg = total degree during follow that was an AM
   AF_deg = total degree during follow that was an AF
   JX_deg = total degree during follow that was a JX
   low_AM_deg = total AM_deg that was low confidence demo assignment
   low_AF_deg = total AF_deg that was low confidence demo assignment
   low_JX_deg = total JX_deg that was low confidence demo assignment
   Follow_length = length of follow in minutes
   Date = date
   Number_SYncs = number of unique synchronized breathing events that occurred
   Group_Size = total individuals present during the focal follow
   Year = year
2. Clean_Degree_Data_SB.csv contains focal follow data from SB. Column descriptions:
   Follow_ID = unique ID for the focal follow conducted
   Dolphin-ID = unique ID of the focal dolphin
   Protocol = type of follow being conducted based on SB definitions
   Demo = Demographic group: AM = Adult Male, AF= Adult Female, AFNN= Mother and neonate calf pair, AFNNN = mother and nonneonate calf pair, JX = juvenile
   Degree = total sync degree during follow
   AM_deg = total degree during follow that was an AM
   AF_deg = total degree during follow that was an AF
   JX_deg = total degree during follow that was a JX
   Follow_length = length of follow in minutes
   Date = date
   Number_SYncs = number of unique synchronized breathing events that occurred
   Group_Size = total individuals present during the focal follow
   Year = year
3. Stranding_Data_background_removed_by_Adult_Sex.csv contains the stranding records during the UME with the background strandings removed for each adult sex class and state based on our methods
4. Stranding_Data_background_removed_by_Age.csv contains the stranding records during the UME with the background strandings removed for each adult sex class and state based on our methods.
5. PC_degree_distributions.csv and SB_degree_distributions.csv are generated using the code "Estimate_degree_distribution_for_infectious_period.R". Each row represents a degree (e.g. row 1 = degree 0, row 2 = degree of 1, etc) and each column (AM, AF, JX) shows the proportion of that demographic that would have that degree over an average DMV infectious period.
6. Mixing_Matrix_data_PC.csv and  Mixing_Matrix_data_SB.csv are generated using the code "Generate_Mixing_Matrices.R". Column descriptions:
   Comb: demo combination for focal:contact
   Mean: average degree for that combination for a focal follow
   Lower: the lower degree estimate of that combination
   Upper: the upper estimate of that combination
   Mean_Ratio: The proportion of all observed syncs that are of that combination
   SD-Ratio: the std of that ratio
   empirical_lower_ratio/upper_ratio: the upper and lower bounds of that proportion

Code Files:


