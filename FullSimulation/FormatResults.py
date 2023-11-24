# defining the libraries
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import csv

M = 500
Scenario = input("Enter Scenario: ")
Topology = input("Enter Topology: (BT/GEANT)")
Pattern = input("Enter Pattern: (random/downstreaming)")
Season = input("Enter Season: (fall/winter/spring/summer)")
#totIntervals = input("Enter Number of Intervals: ")
#totIntervals = int(totIntervals)
# season = input("Enter Season (fall/winter/spring/summer/summerSpecial): ")
if Topology == "BT":
	totNodes = 1008
	totLinks = 3111
else:
	totNodes = 46
	totLinks = 68

Directory = input("Enter the directory of ns-3: (/home/username/) ")

data = []
with open(Directory + "ns-allinone-3.36.1/ns-3.36.1/scratch/CATE/"+ str(Topology) +"Files/Results/DelayValues" + str(Pattern)+ str(Season) + str(Scenario) + ".txt", "r") as file:
	lines_list = file.readlines()
for line in lines_list[0:]:
	data.append([float(val) for val in line.split()][0])

N = len(data) 
# sort the data in ascending order
x = np.sort(data)
# get the cdf values of y
y = np.arange(N) / float(N)

# plotting

#print(len(x), len(y))
Adj = int(len(x)/M)
NewXd = []
NewYd = []
for i in range(M):
	NewXd.append(x[Adj*i])
	NewYd.append(y[Adj*i])
	
#print(NewXd)
#print(NewYd)

plt.plot(x, y, '-k', label = "Sc 1")
plt.plot(NewXd, NewYd, '-g', label = "Sc 1 reduced")

plt.xlabel('Delay (ms)')
plt.ylabel('CDF')
plt.legend()
plt.title('CDF of Delay')
#plt.savefig("DelayCDFMultiple.png")
#plt.show()
plt.figure()

# Get the cdf of hop count values:
data = []
with open(Directory + "ns-allinone-3.36.1/ns-3.36.1/scratch/CATE/"+ str(Topology) +"Files/Results/HopCountValues" + str(Pattern) + str(Season) + str(Scenario) + ".txt", "r") as file:
	lines_list = file.readlines()
for line in lines_list[0:]:
	data.append([float(val) for val in line.split()][0])

N = len(data) 
# sort the data in ascending order
x = np.sort(data)
# get the cdf values of y
y = np.arange(N) / float(N)

# plotting


#print(len(x), len(y))
Adj = int(len(x)/M)
NewXh = []
NewYh = []
for i in range(M):
	NewXh.append(x[Adj*i])
	NewYh.append(y[Adj*i])
	
#print(NewXh)
#print(NewYh)

plt.plot(x, y, '-k', label = "Sc 1")
plt.plot(NewXh, NewYh, '-g', label = "Sc 1 reduced")

plt.xlabel('Hop Count')
plt.ylabel('CDF')
plt.legend()
plt.title('CDF of Hop Count')
#plt.savefig("DelayCDFMultiple.png")
#plt.show()


with open("/home/sawsan/ns-allinone-3.36.1/ns-3.36.1/scratch/CATE/"+ str(Topology) +"Files/Results/CDF" + str(Pattern)+ str(Season) + str(Scenario) + ".csv", 'w', encoding='UTF8') as f:
    writer = csv.writer(f)
    writer.writerow(NewXd)
    writer.writerow(NewYd)
    writer.writerow(NewXh)
    writer.writerow(NewYh)
