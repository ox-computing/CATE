from DijkstraAlgorithm import shortestPath
from HelpingFunctions import GraphTopology, UpdateU, getMaxRate, UpdateCarbon, getUZeroLinks, Graph5, FlowAnalysis, ReadParameters
import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import copy

Topology = input("Enter Topology: (BT/GEANT) ")
Season = input("Enter Season: (fall, winter, spring, summer) ")
directory = input("Enter directory of ns-3: (/home/username/) ")
inputFile1 = directory + "ns-allinone-3.36.1/ns-3.36.1/scratch/CATE/" + Topology + "Files/parameters.txt"
inputFile2 = directory + "ns-allinone-3.36.1/ns-3.36.1/scratch/CATE/" + Topology + "Files/RateMean.txt"

[PortPower, ProcPower, SelectedN, MaxRateST] = ReadParameters(inputFile1, inputFile2)
X = 1
	
PowerSetup = PortPower	
mixed = False
SplittingN = 16 # Used to change the traffic pattern... to be changed by the user...
weightD = 1
weightR = 0.00021
CopyLinks = []

# Get the forecast of carbon intensity per country or region
Forecast = []

with open(directory + "ns-allinone-3.36.1/ns-3.36.1/scratch/"+ Topology +"Files/"+ Season +".txt", "r") as file:
	lines_list = file.readlines()
for line in lines_list[0:]:
	Forecast.append([float(val) for val in line.split()])

#print(Forecast) 
N = len(Forecast)
#Forecast = Forecast[0]
#print("N = ", N)
#print(Forecast)

ECMP = True
FlagNode = []
ratio = 1
Flows1 = []
Flows2 = []
Flows = []
# Generate the graph from the topology file named "filename" with any forecast
graph_BT_original, Clear_G, NodeLevel, NodeLSize, Capacity_G, AllLinks, AllNodes, CoreNodes, MetroNodes, tier_1Nodes, RegionId, AgRate = GraphTopology(directory, Topology, Forecast[0], True, FlagNode, ratio, ProcPower)

# Generate the flows All to All:
for i in AllNodes:
	for j in AllNodes:
		#if (j == i or RegionId[int(i)] == RegionId[int(j)]): continue
		if (j == i): continue
		Flows1.append([i,j, 1])
	
# Generate the flows: mainly downstreaming: every node requires data from the neareset two metro nodes:
for i in CoreNodes:
	neighbors = graph_BT_original[i]
	count = 0
	for edge in neighbors:	
		if edge in MetroNodes:
			Flows2.append([edge, i, 1])
			count += 1
			if count == 2: break
for i in MetroNodes:
	neighbors = graph_BT_original[i]
	count = 0
	for edge in neighbors:	
		if edge in MetroNodes:
			Flows2.append([edge, i, 1])
			count += 1
			break
	if count < 1:
		#print("Not Enough")
		for n in neighbors:
			if count == 1: break
			FurtherNeighbors = graph_BT_original[n]
			for edge in FurtherNeighbors:
				if edge == i: continue
				if edge in MetroNodes:
					Flows2.append([edge, i, 1])
					count += 1
					break
for i in tier_1Nodes:
	neighbors = graph_BT_original[i]
	count = 0
	for edge in neighbors:	
		if edge in MetroNodes:
			Flows2.append([edge, i, 1])
			count += 1
			if count == 2: break
	if count < 2:
		for n in neighbors:
			if count == 2: break
			FurtherNeighbors = graph_BT_original[n]
			for edge in FurtherNeighbors:
				if edge == i: continue
				if edge in MetroNodes:
					Flows2.append([edge, i, 1])
					count += 1
					if count == 2: break		

#print("Number of Flows = ", len(Flows1), len(Flows2))
# We will compare for 3 setups: OSPF, CarbonAware, MST_CarbonAware: 0, 1 and 2.
Carbon_perFlow = [0, 0, 0]
Carbon_LinksEnabled = [0, 0, 0] # this is the embodied carbon emissions caused by enabling links or ports.
CarbonTotal = [0, 0, 0] # this is the sum of embodied and dynamic carbon emissions

FlowLength = [0, 0, 0]
FlowDeltaLength = [0, 0, 0]
NodeFreq = [np.zeros(len(AllNodes)), np.zeros(len(AllNodes)), np.zeros(len(AllNodes))]

ClearG_Untouched = []

MaxRate = [100000000, 100000000, 100000000]
Utilization = [1000, 1000, 1000]
ZeroLinks = [1000, 1000, 1000]
Number_Links_Enb = [0, 0 ,0]
f = open(directory + "ns-allinone-3.36.1/ns-3.36.1/scratch/" + Topology + "Files/PortConf"+ str(PowerSetup) +".txt", "w")

FlowsAvg = 0
for Hour in SelectedN:
	if Hour < SplittingN:
		Flows = Flows1
		FlowsAvg = len(Flows)
		
	elif(not mixed):
		Flows = Flows2
		FlowsAvg = len(Flows)
	else:
		Flows = []
		for i in range(len(Flows2)):
			Flows2[i][2] = weightD
			Flows.append(Flows2[i])
		for i in range(len(Flows1)):
			Flows1[i][2] = weightR
			Flows.append(Flows1[i])
		FlowsAvg = weightR * len(Flows1) + weightD * len(Flows2)
		
	MaxRateS = MaxRateST[Hour]/FlowsAvg
		
	print("Hour: ", Hour)
	for setup in range(3):
		#print("Setup: ", setup)
		if setup == 0: # OSPF no shut down
			graph_BT, Clear_G, NodeLevel, NodeLSize, Capacity_G, AllLinks, AllNodes, CoreNodes, MetroNodes, tier_1Nodes, RegionId, AgRate = GraphTopology(directory, Topology, Forecast[Hour], True, FlagNode, ratio, ProcPower)
			NodeLSize = [int(n)-1 for n in NodeLevel]
		if setup == 1: # Carbon-aware no shut down
			graph_BT, Clear_G, NodeLevel, NodeLSize, Capacity_G, AllLinks, AllNodes, CoreNodes, MetroNodes, tier_1Nodes, RegionId, AgRate = GraphTopology(directory, Topology, Forecast[Hour], False, FlagNode, ratio, ProcPower)
			NodeLSize = [int(n)-1 for n in NodeLevel]
		if setup == 2: # Heuristic
			graph_BT, Clear_G, NodeLevel, NodeLSize, Capacity_G, AllLinks, AllNodes, CoreNodes, MetroNodes, tier_1Nodes, RegionId, AgRate = GraphTopology(directory, Topology, Forecast[Hour], False, FlagNode, ratio, ProcPower)
			NodeLSize = [int(n)-1 for n in NodeLevel]
			g = Graph5(graph_BT, Clear_G, Flows, Forecast[Hour], RegionId, FlagNode, ratio, Capacity_G, NodeLSize, ECMP, MaxRateS, ProcPower, PortPower, X, FlowsAvg)
			graph_BT, DisabledLinks = g.primMST()
			#print(DisabledLinks)
			for link in DisabledLinks:
				try:
					index = AllLinks.index(link)
				except:
					index = AllLinks.index([link[1],link[0]])
				finally:
					f.write(str(index))
					f.write("\t")
			f.write("\n")
		
			
		for S in graph_BT.keys():
			for D in graph_BT[S]:
				Carbon_LinksEnabled[setup] += (Forecast[Hour][RegionId[int(S)]-1] + Forecast[Hour][RegionId[int(D)]-1]) * PortPower / 2 # in Ws * gCO2/KWh so should be /3600000 to become in gCO2 but keep it now to keep the accuracy
		
		length = sum([len(graph_BT[key]) for key in graph_BT.keys()]) / 2
		#print(length)
		Number_Links_Enb[setup] = length
		

		Clear_G, CarbonFlows, AverageFlowLength, NodeFreq1, q = FlowAnalysis(graph_BT, copy.deepcopy(Clear_G), Flows, FlagNode, Forecast[Hour], RegionId, ratio, ECMP, NodeLSize, ProcPower)
		Carbon_perFlow[setup] = CarbonFlows
		FlowLength[setup] = AverageFlowLength
		for index in range(len(NodeFreq1)):
			NodeFreq[setup][index] += NodeFreq1[index] / N
			
		max_rate = getMaxRate(Clear_G, Capacity_G)
		MaxRate[setup] = max_rate
		#if max_rate < MaxRate[setup]: MaxRate[setup] = max_rate
		U, ZL = getUZeroLinks(graph_BT, Clear_G, max_rate, Capacity_G)
		if U < Utilization[setup]: Utilization[setup] = U
		if ZL < ZeroLinks[setup]: ZeroLinks[setup] = ZL
		#print(MaxRate[setup], MaxRateS)
		#print("Rate = ", MaxRate[setup], " and max Throughput = ", MaxRate[setup]*FlowsAvg)
		#print("Actual Throughput = ", MaxRateS*FlowsAvg, "Gbps, U =", U, "%", "Carbon of Ports =", Carbon_LinksEnabled[setup]/3600000, "gCO2, Carbon Dynamic =", Carbon_perFlow[setup] * MaxRateS * 1000/3600000, "gCO2")
		#print("Net Improvement w.r.t OSPF =", (Carbon_LinksEnabled[0]+Carbon_perFlow[0] * MaxRateS * 1000 - Carbon_LinksEnabled[setup] - Carbon_perFlow[setup] * MaxRateS * 1000)/3600000 * 30*60, "gCO2 per 30 min, Percentage Improvement=", (Carbon_LinksEnabled[0]+Carbon_perFlow[0] * MaxRateS * 1000 - (Carbon_LinksEnabled[setup]+Carbon_perFlow[setup] * MaxRateS * 1000))/(Carbon_LinksEnabled[0]+Carbon_perFlow[0] * MaxRateS * 1000) * 100, "%, Improvement in Dynamic =", (Carbon_perFlow[0] - Carbon_perFlow[setup])/(Carbon_perFlow[0]) * 100, "%, Improvement in Port-Prop =", (Carbon_LinksEnabled[0] - Carbon_LinksEnabled[setup])/(Carbon_LinksEnabled[0]) * 100, "%, Number of Links Enabled =", Number_Links_Enb[setup]/Number_Links_Enb[0]*100, "%", "Number of Links Disabled =", (Number_Links_Enb[0]-Number_Links_Enb[setup])/Number_Links_Enb[0]*100, "%")
		

		
		
		
		
		
		
		
		


