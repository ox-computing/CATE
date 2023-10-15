import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
from matplotlib.colors import Normalize
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


Directory = input("Enter the directory of ns-3: (/home/username/) ")
Topology = input("Enter the Topology: (BT/GEANT) ")
Pattern = input("Enter the pattern: (random/downstreaming/mixed) ")
Season = input("Enter Season: (fall/winter/spring/summer) ")
Scen = [int(input("Enter Scenario to be included: (1-9,20) "))]
Resp = input("Do you want to include another Scenario? (y/n) ")
while(Resp == "y"):
	Scen.append(int(input("Enter Scenario to be included: (1-9,20) ")))
	Resp = input("Do you want to include another Scenario? (y/n) ")	

LABELS = ['OSPF', 'IncD', 'C', 'C+IncD', 'C+Ptyp', 'C+E-label', 'Ptyp', 'E-label', 'CE', 'CATE']


Scenarios = len(Scen)

CFixed = []
CPorts = []
m1 = 0
m2 = 0
scale = 1


intervals = 9


DATA = []
for s in Scen:  
	data = []
	with open(Directory + "ns-allinone-3.36.1/ns-3.36.1/scratch/CATE/"+Topology+"Files/Results/Results" + str(Pattern)+ str(Season) + str(s) + ".txt", "r") as file:
		lines_list = file.readlines()
	#print(lines_list)
	for line in lines_list:
		data.append([float(val) for val in line.split()])
	DATA.append(data)

for s in range(len(Scen)):
	if Scen[s] == 20:
		Scen[s] = 9;
	else:
		Scen[s] = Scen[s] -1
		
# Draw carbon savings % 

OverallImpC = []
for s in range(Scenarios):
	v = []
	for i in range(intervals):
		c = DATA[s][2][i] + DATA[s][7][i] + DATA[s][9][i]
		c0 = DATA[0][2][i] + DATA[0][7][i] + DATA[0][9][i]
		v.append((c0-c)/c0 * 100)
	OverallImpC.append(v)
	
plt.rcParams['font.size'] = 14
n = 2
w = 10
for s in range(len(Scen)):
	if s > 0:
		plt.plot(range(intervals), OverallImpC[s],  label=LABELS[Scen[s]], linewidth=n, markersize=w)
plt.xticks([0, 2, 4, 6, 8], ['9:00', '11:00', '13:00', '15:00', '17:00']) # Adjust as needed
plt.legend()
plt.grid()
plt.xlim(0, intervals-1)
plt.xlabel("Time")
plt.ylabel("% Carbon Savings w.r.t OSPF")
plt.show()

OverallImpE = []
for s in range(Scenarios):
	v = []
	for i in range(intervals):
		c = DATA[s][1][i] + DATA[s][6][i] + DATA[s][8][i]
		c0 = DATA[0][1][i] + DATA[0][6][i] + DATA[0][8][i]
		v.append((c0-c)/c0 * 100)
	OverallImpE.append(v)
	

for s in range(len(Scen)):
	if s > 0:
		plt.plot(range(intervals), OverallImpE[s],  label=LABELS[Scen[s]], linewidth=n, markersize=w)
plt.xticks([0, 2, 4, 6, 8], ['9:00', '11:00', '13:00', '15:00', '17:00'])
plt.legend()
plt.grid()
plt.xlim(0, intervals-1)
plt.xlabel("Time")
plt.ylabel("% Energy Savings w.r.t OSPF")
plt.show()


# Draw Heatmap for number of flows per region (BT)
Y = []
max1 = 0
for s in range(Scenarios): 
	arr = np.array(DATA[s][5])
	Y.append(arr)
	

Ylabel = []
for s in range(len(Scen)):
	Ylabel.append(LABELS[Scen[s]])
	
# Plot the heatmap
plt.rcParams['font.size'] = 10.5
plt.figure(figsize=(16,3))
#cmap = sns.diverging_palette(130, 10, sep=1, as_cmap=True)
cmap = LinearSegmentedColormap.from_list('red-blue', ['lightblue', 'red'])
if Topology == "GEANT":
	Xlabel= ["Iceland", "Ireland", "UK", "Portugal", "Spain", "France", "Netherlands", "Belgium", "Luxembourg", "Switzerland", "Italy", "Norway", "Denmark", "Germany", "Finland", "Sweden", "Estonia", "Latvia", "Lithiuania", "Poland", "Czechia", "Ukraine", "Slovakia", "Austria", "Slovenia", "Croatia", "Hungary", "Moldova", "Romania", "Serbia", "Albania", "Montenegro", "Bulgaria", "North Macedonia", "Turkey", "Greece", "Malta", "Cyprus", "Israel"]
	Xlabel2= ["IS", "IE", "UK", "PT", "ES", "FR", "NL", "BE", "LU", "CH", "IT", "NO", "DK", "DE", "FI", "SE", "EE", "LV", "LT", "PL", "CZ", "UA", "SK", "AT", "SI", "HR", "HU", "MD", "RO", "RS", "AL", "ME", "BG", "MK", "TR", "EL", "MT", "CY", "IL"]
	heat_map = sns.heatmap( Y, cmap = cmap, linewidth = 5 , annot = True, annot_kws={'rotation': 90}, fmt='.3g', xticklabels=Xlabel2, yticklabels=Ylabel)
else:
	heat_map = sns.heatmap( Y, cmap = cmap, linewidth = 5 , annot = True, annot_kws={'rotation': 90}, fmt='.3g', yticklabels=Ylabel)


plt.title( "Average Intensity of Flows per Region (Gbps)" )
plt.xticks(rotation=90)
plt.show()

