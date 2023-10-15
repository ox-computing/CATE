import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Plot the overall rate for BT topology (Fig. 8)
x = x = np.linspace(0, 48, 49)
y = [0.41, 0.35, 0.3, 0.25, 0.2, 0.18, 0.16, 0.155, 0.15, 0.19, 0.26, 0.34, 0.42, 0.4, 0.39, 0.39, 0.38, 0.38, 0.39, 0.4, 0.4, 0.4, 0.4, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.47, 0.485, 0.5, 0.54, 0.59, 0.62, 0.7, 0.74, 0.82, 0.88, 0.92, 0.96, 1, 0.98, 0.9, 0.8, 0.7, 0.6, 0.5, 0.41]

plt.rcParams['font.size'] = 14
plt.plot(x, y, color='g', alpha=.5)
plt.xlim(0, 48)
plt.ylim(0, 1)
plt.xlabel("Time (h)")
plt.ylabel("Normalized Rate")
plt.xticks([0, 12, 24, 36, 48], ['12AM', '6AM', '12PM', '6PM', '12AM'])
plt.fill_between(x, y, 0, color='g', alpha=.5)
plt.savefig("TrafficPattern.png")
plt.show()

# Plot the National Carbon Intensity in the UK (Fig. 1)
w = [86, 77, 76, 75, 81, 74, 79, 73, 74, 74, 73, 80, 83, 102, 110, 133, 140, 153, 156, 146, 139, 137, 134, 119, 113, 107, 108, 107, 114, 114, 119, 123, 142, 141, 142, 141, 143, 146, 149, 148, 143, 147, 152, 144, 128, 111, 104, 91, 90]

sp = [252, 252, 252, 249, 245, 238, 235, 242, 239, 244, 245, 246, 252, 253, 244, 234, 226, 219, 219, 218, 219, 209, 202, 198, 199, 200, 201, 186, 193, 203, 197, 199, 204, 209, 215, 221, 218, 222, 229, 229, 237, 232, 231, 231, 241, 233, 234, 235, 238]

sm = [226, 230, 230, 238, 231, 228, 231, 225, 218, 218, 232, 234, 233, 235, 239, 222, 221, 213, 212, 204, 201, 199, 193, 192, 188, 185, 180, 173, 175, 173, 174, 181, 184, 191, 192, 192, 194, 199, 201, 203, 207, 213, 209, 211, 208, 207, 204, 199, 194]

f = [84, 69, 76, 77, 74, 76, 73, 74, 78, 79, 80, 73, 82, 94, 101, 114, 111, 98, 86, 79, 76, 87, 81, 82, 84, 85, 86, 87, 96, 110, 119, 124, 142, 153, 160, 162, 164, 146, 142, 144, 150, 142, 122, 114, 97, 80, 80, 68, 67]

x = np.linspace(0, 48, 49)

plt.plot(x, w, label = 'Winter')
plt.plot(x, sp, label = 'Spring')
plt.plot(x, sm, label = 'Summer')
plt.plot(x, f, label = 'Fall')

plt.xlim(0, 48)
plt.xticks([0, 12, 24, 36, 48], ['12AM', '6AM', '12PM', '6PM', '12AM'])
plt.xlabel("Time (h)")
plt.ylabel("UK National Carbon Intensity (gCO2/kWh)")
plt.legend()
plt.savefig("CIPattern.png")
plt.show()


# Plot the Regional Carbon Intensity in the UK (Fig. 2)
Scotland = []
NorthWales = []
London = []

directory = input("Enter the ns-3 directory: (/home/username/) ")

with open(directory + "ns-allinone-3.36.1/ns-3.36.1/scratch/project-carbon-aware-routing-sim/BTFiles/fall (copy).txt", "r") as file:
	lines_list = file.readlines()
	for line in lines_list[0:]:
		Scotland.append([float(val) for val in line.split()][0])
		NorthWales.append([float(val) for val in line.split()][5])
		London.append([float(val) for val in line.split()][12])
		
x = np.linspace(0, 47, 48)

plt.plot(x, Scotland, label = 'Scotland')
plt.plot(x, NorthWales, label = 'NorthWales')
plt.plot(x, London, label = 'London')


plt.xlim(0, 47)
plt.xticks([0, 12, 24, 36], ['12AM', '6AM', '12PM', '6PM'])
plt.xlabel("Time (h)")
plt.ylabel("UK Regional Carbon Intensity (gCO2/kWh)")
plt.legend()
plt.savefig("CIPattern.png")
plt.show()

