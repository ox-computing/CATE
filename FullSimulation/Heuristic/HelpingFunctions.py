import sys
import numpy as np
import copy
from DijkstraAlgorithm import shortestPath
from multiprocessing import Process, Queue, Pool

# Trandforming Inet Topology to a graph
def GraphTopology(directory, filename, Forecast, OSPF, FlagNode, ratio, ProcPower):
	G = {} # Graph of the topology
	Clear_G = {}
	Capacity_G = {}
	data = []
	RegionId = []
	NodeLevel = []
	linkRate = []
	Edges = []
	tier_1Nodes = []
	CoreNodes = []
	MetroNodes = []
	AllNodes = []
	AgRate = []
	AllLinks = []
	NodeLSize = []
	
	with open(directory + "ns-allinone-3.36.1/ns-3.36.1/scratch/CATE/"+ filename +"Files/Topology.txt", "r") as file:
		lines_list = file.readlines()
	
	for line in lines_list[0:]:
		data.append([val for val in line.split()])
	
	totNodes = int(data[0][0])
	totLinks = int(data[0][1])
	
	Links_perNode = np.zeros(totNodes)
	
	#print(totNodes, totLinks)
	for i in range(1, 1+ totNodes):
		G[data[i][0]] = {}
		Clear_G[data[i][0]] = {}
		Capacity_G[data[i][0]] = {}
		RegionId.append(int(data[i][1]))
		NodeLevel.append(data[i][2])
		AgRate.append(0)
		AllNodes.append(data[i][0])
		if (data[i][2] == '3'): tier_1Nodes.append(data[i][0])
		if (data[i][2] == '2'): MetroNodes.append(data[i][0])
		if (data[i][2] == '1'): CoreNodes.append(data[i][0])
	
	for i in range(1+totNodes, 1+totNodes+totLinks):
		From = data[i][0]
		To = data[i][1]
		AllLinks.append([From, To])
		if (int(From) in FlagNode):
			costFrom = 1 + Forecast[RegionId[int(From)]-1] * ratio * ProcPower[int(NodeLevel[int(From)])-1] * 1000
		else: 
			costFrom = 1 + Forecast[RegionId[int(From)]-1] * ProcPower[int(NodeLevel[int(From)])-1] * 1000
		if (int(To) in FlagNode):
			costTo = 1 + Forecast[RegionId[int(To)]-1] * ratio * ProcPower[int(NodeLevel[int(From)])-1] * 1000
		else:
			costTo = 1 + Forecast[RegionId[int(To)]-1] * ProcPower[int(NodeLevel[int(From)])-1] * 1000
		if OSPF:
			costFrom = 1
			costTo = 1
		G[From][To] = costTo
		G[To][From] = costFrom
		Clear_G[From][To] = 0
		Clear_G[To][From] = 0
		Capacity_G[From][To] = data[i][2]
		Capacity_G[To][From] = data[i][2]
		AgRate[int(From)] += float(data[i][2])
		AgRate[int(To)] += float(data[i][2])
		Links_perNode[int(From)] += 1
		Links_perNode[int(To)] += 1
		#linkRate.append(data[i][2])
		#Edges.append([From, To])
	for i in range(totNodes):
		if filename == "GEANT":
			NodeLSize.append(3)
		elif Links_perNode[i] <= 32:
			NodeLSize.append(2)
		elif Links_perNode[i] <= 48:
			NodeLSize.append(1)
		else:
			NodeLSize.append(0)
		
		
	#print(G)
	return G, Clear_G, NodeLevel, NodeLSize, Capacity_G, AllLinks, AllNodes, CoreNodes, MetroNodes, tier_1Nodes, RegionId, AgRate


def UpdateU(Clear_G, list_P, weight):
	value = weight / len(list_P)
	for P in list_P:	
		for i in range(len(P)-1):
			Clear_G[P[i]][P[i+1]] += value
	return Clear_G

def getMaxRate(Clear_G, Capacity_G):
	r = 10000
	for node in Clear_G:
		for edge in Clear_G[node]:
			if float(Clear_G[node][edge]) != 0:
				value = float(Capacity_G[node][edge])/float(Clear_G[node][edge])
				if value < r:
					r = value
	return r
	
	
def UpdateCarbon(NodeLoad, list_P, RegionId, Forecast, AgRate, rate, weight):
	carbon = 0
	SumValue = 0
	value = rate * weight / len(list_P)
	for P in list_P:
		for node in P:
			carbon += value * Forecast[RegionId[int(node)]-1] / AgRate[int(node)]
			NodeLoad[int(node)] += 1
	return carbon, NodeLoad
	
def getUZeroLinks(G, Clear_G, rate, Capacity_G):
	tot = 0
	num = 0
	load = 0
	MaxCap = 0
	for node in G:
		for edge in G[node]:
			tot += 1
			load += float(Clear_G[node][edge]) * rate
			MaxCap += float(Capacity_G[node][edge])
			if Clear_G[node][edge] == 0: num += 1
	U = 100 * load / MaxCap
	return U, num/tot

 
def run_cpu_tasks_in_parallel(Flows, graph, ECMP):
    inputs = []
    for i in range(len(Flows)):
    	inputs.append([Flows[i], graph, ECMP])
    with Pool() as p:
    	q = p.map(ProcessFlow, inputs)
    #print("q is completed!")
    return q
    
def ProcessFlow(inputs):
	#print(inputs[0])
	list_P = shortestPath(inputs[2], inputs[1], inputs[0][0], inputs[0][1])
	return [list_P, inputs[0][2]]

def FlowAnalysis(graph, ClearG, Flows, FlagNode, Forecast, RegionId, ratio, ECMP, NodeLSize, ProcPower):
	FlowLength = []
	Carbon_perFlow = []
	NodeFreq = np.zeros(len(graph))
	q = run_cpu_tasks_in_parallel(Flows, graph, ECMP)
	#counter = 0
	for list_P1 in q:
		list_P = list_P1[0]
		weight = list_P1[1]
		ClearG = UpdateU(ClearG, list_P, weight)
		carbon = 0
		FlowLength.append((len(list_P[0])-1))
		for P in list_P:
			for c in P:
				if (int(c) in FlagNode):
					carbon += Forecast[RegionId[int(c)]-1]/len(list_P) * ratio * ProcPower[int(NodeLSize[int(c)])]
				else:
					carbon += Forecast[RegionId[int(c)]-1]/len(list_P) * ProcPower[int(NodeLSize[int(c)])]
				NodeFreq[int(c)] += weight/len(list_P)
		Carbon_perFlow.append(weight * carbon)
		#counter += 1
		#print("Completed:", counter/len(q) * 100, "%")
		
	CarbonFlows = sum(Carbon_perFlow)
	AverageFlowLength = sum(FlowLength)
	return ClearG, CarbonFlows, AverageFlowLength, NodeFreq, q

def ReadParameters(inputFile1, inputFile2):
	data = []
	with open(inputFile1, "r") as file:
		lines_list = file.readlines()
	for line in lines_list[0:]:
		data.append([val for val in line.split()])
		
	N = int(data[3][2])
	L = int(data[1][2])
	PortPower = float(data[9][2])
	ProcPower = []
	for i in range(L):
		ProcPower.append(float(data[13][2+i]))
	SelectedN = range(int(data[3][2]))

	with open(inputFile2, "r") as file:
		lines_list = file.readlines()
	line = lines_list[0:][0]
	Rate = np.array([float(val) for val in line.split()])
	MaxRateST = Rate * float(data[15][2])
	return [PortPower, ProcPower, SelectedN, MaxRateST] 
		
class Graph5():
	def __init__(self, G, ClearG, Flows, Forecast, RegionId, FlagNode, ratio, Capacity_G, NodeLSize, ECMP, Rate0, ProcPower, PortPower, X, FlowsAvg):
		self.V = len(G)
		self.ClearG = ClearG
		self.CarbonGraph = G # carbon intensity information 
		self.Flows = Flows
		self.Forecast = Forecast
		self.FlagNode = FlagNode
		self.ratio = ratio
		self.RegionId = RegionId
		self.Capacity_G = Capacity_G
		self.NodeLSize = NodeLSize
		self.ECMP = ECMP
		self.FlowPaths = [[] for _ in range(len(self.Flows))]
		self.Rate0 = Rate0
		self.ProcPower = ProcPower
		self.PortPower = PortPower
		self.X = X
		self.FlowsAvg = FlowsAvg
		
	def AdjustGraphFormat(self, graph):
		newG = {}
		for i in range(self.V):
			newG[str(i)] = {}
			
		for key1 in graph.keys():
			for key2 in graph[key1].keys():
				#print(key1, key2)
				#print(self.CarbonGraph[key1][key2], self.CarbonGraph[key2][key1], self.FlowGraph[key1][key2],  self.FlowGraph[key2][key1])
				if (self.FlowGraph[key1][key2] + self.FlowGraph[key2][key1] ==0):
					newG[key1][key2] = sys.maxsize
				else: 
					newG[key1][key2]= 1 + 1000 * (self.CarbonGraph[key1][key2]+ self.CarbonGraph[key2][key1]) / (self.FlowGraph[key1][key2] + self.FlowGraph[key2][key1])
	
		return newG

	def CarbonLinks(self, graph):
		Carbon_LinksEnabled = 0
		for S in graph.keys():
			for D in graph[S]:
				Carbon_LinksEnabled += (self.Forecast[self.RegionId[int(S)]-1] + self.Forecast[self.RegionId[int(D)]-1]) * self.PortPower / 2 # divide by 2 because every value is added twice
		return Carbon_LinksEnabled

	def FormatResult(self, graph):
		G = {} # Graph of the topology
		for i in range(self.V):
			G[str(i)] = {}
		
		for key1 in graph.keys():
			for key2 in graph[key1].keys():
				G[key1][key2] = self.CarbonGraph[key1][key2]
		return G
		
	def GraphStats(self, G, FlowGraph):
		max_rate = getMaxRate(FlowGraph, self.Capacity_G)
		U, ZL = getUZeroLinks(G, FlowGraph, max_rate, self.Capacity_G)
		return max_rate, U, ZL
		
	def primMST(self):
		DisabledLinks = []
		ClearG, CarbonFlows, AverageFlowLength, NodeFreq, FlowPaths1= FlowAnalysis(self.CarbonGraph, copy.deepcopy(self.ClearG), self.Flows, [], self.Forecast, self.RegionId, self.ratio, True, self.NodeLSize, self.ProcPower)
		#print("Start: ", self.CarbonGraph)
		self.FlowGraph = ClearG
		self.FlowPaths = FlowPaths1
		#print("Initial flows computed")
		#print("Flows before: ", self.FlowGraph)
		self.graph = self.AdjustGraphFormat(self.CarbonGraph)
		self.baselineFlows = CarbonFlows
		self.baselineLinks = self.CarbonLinks(self.CarbonGraph)
		SetImprovement = 0
		NecessaryLinks = []
		more = True
		while(True):
			if not more: break
			newG = copy.deepcopy(self.graph)
			length = sum([len(newG[key]) for key in newG.keys()]) / 2
			#print(newG)
			links = []
			
			#for count in range(max(1, int((length-self.V)/3))):
			for count in range(self.X):
				max1 = -1
				for S in newG.keys():
					for D in newG[S]:
						value = newG[S][D]
						#print("value", value)
						if value > max1 and len(newG[S])>1 and len(newG[D])>1 and ([S,D] not in NecessaryLinks) and ([D,S] not in NecessaryLinks): 
							max1 = value
							link = [S, D]
							#print(link)
							#links.append(link)
				
				if max1 == -1: 
					#print("No more links can be shut down!")
					more = False
					break
				#print(link, value)
				del newG[link[0]][link[1]]
				del newG[link[1]][link[0]]
				links.append(link)
			if len(links) == 0:
				break
			newG = self.FormatResult(newG)
			
			try: 
				#print("Number of Links selected = ", len(links))
				ClearG, CarbonFlows, AverageFlowLength, NodeFreq, FlowPaths1 = FlowAnalysis(newG, copy.deepcopy(self.ClearG), self.Flows, [], self.Forecast, self.RegionId, self.ratio, True, self.NodeLSize, self.ProcPower)

			except: 
				#print("Graph not connected!")
				for link in links:
					NecessaryLinks.append(link)
				#print("Necessary Links: ", len(NecessaryLinks))
				#print(NecessaryLinks)
			else:
				#ImprovementLinks = (self.baselineLinks - self.CarbonLinks(newG))/self.baselineLinks * 100
				#ImprovementFlows = (self.baselineFlows - CarbonFlows)/ self.baselineFlows * 100
				#Improvement = ImprovementLinks * self.FixRatio + ImprovementFlows * (1 - self.FixRatio)
				#print(self.baselineLinks, self.CarbonLinks(newG), self.baselineFlows*self.Rate0*1000, CarbonFlows * self.Rate0 * 1000)
				#print(self.baselineLinks - self.CarbonLinks(newG))
				#print(self.baselineFlows*self.Rate0*1000 - CarbonFlows * self.Rate0 * 1000)
				Improvement = (self.baselineLinks+self.baselineFlows*self.Rate0*1000)-(self.CarbonLinks(newG) + CarbonFlows * self.Rate0 * 1000)
				#print("Improvement = ", Improvement/3600000 * 30*60, "gCO2 per 30 min")
				max_rate, U, ZL = self.GraphStats(newG, ClearG)
				#print("New max rate = ", max_rate, self.Rate0, max_rate* self.FlowsAvg)
				
				if Improvement > SetImprovement and max_rate >= self.Rate0:
					self.FlowGraph = ClearG
					self.FlowPaths = FlowPaths1
					#print("Flows After: ", self.FlowGraph)
					self.graph = self.AdjustGraphFormat(newG)
					SetImprovement = Improvement
					#print("Improvement is higher and graph is connected")
					for link in links:
						DisabledLinks.append(link)
				else:
					#print("Improvement dimished! OR max rate unreachable! BREAK")
					break
		
		return self.FormatResult(self.graph), DisabledLinks
