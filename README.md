# CATE
## Introduction
This is a project repository for carbon aware routing simulations using the ns-3 simulator tool. 
Carbon-aware routing is introducing the knowledge of the carbon intensity into the routing process. It aims at reducing the overall carbon emissions associated with routers in a network. Intra-domain routing is considered in this project which is routing within one Autonomous System (AS).
The full description of the simulation setup, metrics and topologies can be found in the paper.
This simulation supports the evaluation of the benefits of carbon-aware routing based on two approaches: (1) Changing the link costs based on the metrics that are defined in the paper, (2) Changing the link costs and shutting down links to further reduce the overall carbon emissions. 

## Installation
The following steps are needed to run this code:

1. Install the ns-3 simulator. The version used in this project is 3.36.1. The following YouTube video is an easy tutorial for installing ns3 in Ubuntu: https://www.youtube.com/watch?v=3lWeCGPiWWM. Configure ns-3 into optimized mode.
2. Clone this repository under the directory ~/ns-allinone-3.36.1/ns-3.36.1/scratch.

## Structure
The repository is strutured as follows:
1. "FullSimulation" folder: this is the main folder containing all the codes. These codes require a topology and a set of parameters.
2. "GEANTFiles" folder: an example folder containing the required information about the GEANT topology.
3. "BTFiles" folder: Similar structure to "GEANTFiles" but related to the BT topology. Unfortunately, the BT topology cannot be made public but we include the simulation parameters for reference.

The "FullSimulation" folder includes:
1. "NoShutDown" folder: contains "FullSimulation.cc" script that changes the link costs based on the metrics defined in the paper.
2. "LinkShutDown" folder: contains "FullSimulationP.cc" that further shuts down links.
3. "Heuristic" folder: contains the theoretical analysis for CATE. The main file under this directory is "HeuristicTest.py" while other scripts are helping functions.
4. "FormatResults.py" script: used to format the results after running a simulation.
5. "CarbonEnergyFigures.py" script: used to plot the figures for the resulting carbon and energy savings for every scenario.
6. "CDFFigures.py" script: used to plot the CDF figures for the delay and the hop count for every scenario. 

Every topology file (i.e. "GEANTFiles") includes:
1. A topology file in the INET format: "Topology.txt".
2. The normalized average rate per flow for every interval of time: "RateMean.txt"
3. The carbon intensity values pr region for every interval and for every season: "fall.txt", "winter.txt", "spring.txt" and "summer.txt".
4. A script for the simulation parameters to be read by the ns-3 code: "parameters.txt".
These files should be provided for every topology to be evaluated. The topology folder should be named as "TOPOLOGYFiles", for example for Geant the folder is named "GEANTFiles".

## Usage
1. Fill the parameters in the topology file.
2. Run the first approach while varying the scenario, the mapping of scenarios is as follows: 1 = OSPF, 2 = Incremental dynamic power per unit traffic only, 3 = Carbon Intensity only, 4 = Carbon Intensity + Incremental dynamic power per unit traffic, 5 = Typical Power + Carbon Intensity, 6 = Energy Labelling + Carbon Intensity, 7 = Typical Power only, 8 = Energy Labelling only, 9 = Carbon Emissions. You need to specify the topology, the season and the scenario when running the code. Moreover, you need to specify the directory at which you installed ns-3, mainly "/home/username/". You can change it in the code or give it as an input argument. The command to run is for example: ./ns3 run "scratch/CATE/FullSimulation/NoShutDown/FullSimulation.cc --Topology=GEANT --Season=winter --Scenario=4 --Directory=/home/username/"
3. Run the Heuristic script to get the number of links to shut down per interval of time. You will be prompted to enter the topology name. The command to run is: python3 scratch/CATE/FullSimulation/Heuristic/HeuristicTest.py and then insert the topology name all in uppercase, example: "GEANT".
4. Run the second approach that corresponds to CATE. The scenario number for this approach is 20 by default. The scenario number if used later to facilitate plotting graphs. You need to specify the topology and the season in the command line as well as the directory. The command to run is for example: ./ns3 run "scratch/CATE/FullSimulation/LinkShutDown/FullSimulationP.cc --Topology=GEANT --Season=winter --Directory=/home/username/"

## License


## Citation
When referencing this work, please use the following citation:

S. El Zahr, P. Gunning, and N. Zilberman. 2023. "Exploring the Benefits of Carbon-Aware Routing". In Proceedings of the 19th International Conference on emerging Networking EXperiments and Technologies (CoNEXT '23).
