import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
from matplotlib.colors import Normalize
import seaborn as sns
import csv

n = 2
Directory = input("Enter the directory of ns-3: (/home/username/)")
Topology = input("Enter the Topology: (BT/GEANT) ")
Pattern1 = input("Enter the pattern: (random/downstreaming/mixed) ")
Season = input("Enter Season: (fall/winter/spring/summer) ")
Scenario = [int(input("Enter Scenario to be included: (1-9,20) "))]
Resp = input("Do you want to include another Scenario? (y/n)")
while(Resp == "y"):
	Scenario.append(int(input("Enter Scenario to be included: (1-9,20) ")))
	Resp = input("Do you want to include another Scenario? (y/n)")	

LABELS = ['OSPF', 'IncD', 'C', 'C+IncD', 'C+Ptyp', 'C+E-label', 'Ptyp', 'E-label', 'CE', 'CATE']
DATA1 = []

for s in Scenario:
	print(s)
	data = []
	with open(Directory + 'ns-allinone-3.36.1/ns-3.36.1/scratch/CATE/'+ Topology +'Files/Results/CDF'+ str(Pattern1)+ str(Season) + str(s)+ '.csv', newline='') as csvfile:
		spamreader = csv.reader(csvfile, quotechar='|')
		for row in spamreader:
			data.append(row)  
	data = np.array(data)
	data = data.astype(np.float64)
	DATA1.append(data)

for s in range(len(Scenario)):
	if Scenario[s] == 20:
		Scenario[s] = 9;
	else:
		Scenario[s] = Scenario[s] -1
		
# Delay Curves	  
plt.rcParams['font.size'] = 14
w = 10
n = 2
for s in range(len(Scenario)):
	plt.plot(DATA1[s][0], DATA1[s][1], linewidth=n, label=LABELS[Scenario[s]], markersize=w)

plt.xlabel('Delay (ms)')
plt.ylabel('CDF')
plt.xlim(0, 150)
plt.ylim(0, 1)
plt.legend()
plt.grid()
plt.show()


# Hop Count Curves
for s in range(len(Scenario)):
	plt.plot(DATA1[s][2], DATA1[s][3], linewidth=n, label=LABELS[Scenario[s]], markersize=w)

plt.xlabel('Hop Count')
plt.ylabel('CDF')
plt.xlim(1, 11)
plt.ylim(0, 1)
plt.legend()
plt.grid()
plt.show()
