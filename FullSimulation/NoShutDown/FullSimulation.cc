#include <ctime>
#include <sstream>
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/netanim-module.h"
#include "ns3/topology-read-module.h"
#include <list>
#include "ns3/flow-monitor-module.h"
#include <numeric>
#include "ns3/traffic-control-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("CarbonAwareProject");
std::string dir = "/home/username/";
std::string Topology = "GEANT";
std::string Pattern = "random";
// The following variables are initiated here but the values are changed later in the code.
double Sec_per_Interval = 0.1; // number of seconds to consider per interval
int Number_of_Intervals = 9; // NUmber of time intervals (usually based on carbon intensity intervals)
double Interval_Adjustment = Sec_per_Interval; // to map every second of simulation to the interval
double Rscale = 1; // For reduced simulation time, if needed
double rate_max = 35/Rscale; //this rate is per flow and will be multipied by a normalized mean at every interval
int nb_Applications = 10; // Number of parallel flows between a pair of nodes // Multiple smaller flows rather than a big flow that is not transmitted because of RTT limitations in TCP
double SimStart = 0;
int CustomPacketSize = 1500; // 64 KB = 65535 bits is the maximum value, can be increased to reduce the simulation time.
double offset = 0.95; // offset time in seconds to let the system set up before starting the transmission of flows.
double PortPower = 4.5; // in W

std::vector<double> DPRating = {15}; // This is 10000 times lamda used only for scaling
std::vector<double> ERating = {10}; // The corresponding energy rating according to the description in the paper.
std::vector<double> TypPower = {187}; // Typical Power in W
std::vector<double> lamda = {0.0015}; // this is the power consumption per unit of traffic in Watt/Mbps for every node level
std::vector<double> FixedNPower = {41}; // Static Power in W

int l0 = 0;
NodeContainer nodes; 
int totRegions;
std::vector<int> RegionIds; 
std::vector<double> HeatR;  // Heat map per region that shows the number of flows in every region (not every node)
std::vector<double> AgRate; // in Gbps
std::vector<int> NodeLevel;
std::vector<int> CoreNodes;
std::vector<int> MetroNodes;
std::vector<int> Tier1Nodes;
std::vector<int> AllNodes;
std::vector<std::vector<int>> Neighbors;
std::vector<std::vector<int>> Flows;
std::vector<NodeContainer> nc; 
std::vector<std::vector<int>> LinkRegions;
int totNodes;
int totLinks;
std::vector<double> weight; // in Gbps
double Gscale = 1000000000; // to adjust the values of Gbps // do not touch
std::vector<std::string> delay_perLink;  // in ns

double Throughput = 0;
double ThroughputR = 0; // for random traffic
double ThroughputD = 0; // for downstreaming traffic
std::vector<double> Throughput_per_Interval;
std::vector<double> load; // Number of packets passing through a node
std::vector<double> load_prev; 
std::vector<double> load_per_int;
int dropped_packets = 0;
std::vector<double> linkLoad;
std::vector<double> Max_linkLoad;
std::vector<double> NodeLoad;
std::vector<Ipv4InterfaceContainer> ipic;
std::vector<NetDeviceContainer> ndc;

double EnergyConsumption = 0;
double CarbonEmissions = 0;
std::vector<double> EnergyConsumption_per_intervalD;
std::vector<double> CarbonEmissions_per_intervalD;
std::vector<double> EnergyConsumption_per_intervalP;
std::vector<double> CarbonEmissions_per_intervalP;
std::vector<double> EnergyConsumption_per_intervalF;
std::vector<double> CarbonEmissions_per_intervalF;
std::vector<double> Embedded_EC_perNode;
std::vector<double> Operational_EC_perNode;
std::vector<double> Embedded_CE_perNode;
std::vector<double> Operational_CE_perNode;

std::vector<double> hopCount_array;
std::vector<double> delay_array;

int nbAppRunning = 0;
std::vector<int> Possible_size {52, 56, 60, 68, 76, 1500, 6000, 65535};
std::vector<int> Nb_Gen_perSize {0, 0, 0, 0, 0, 0, 0, 0};

// Special class of tag that allows to record the time when a packet was sent by the socket.
class TimestampTag : public Tag {
public:
  static TypeId GetTypeId (void);
  virtual TypeId GetInstanceTypeId (void) const;

  virtual uint32_t GetSerializedSize (void) const;
  virtual void Serialize (TagBuffer i) const;
  virtual void Deserialize (TagBuffer i);

  // these are our accessors to our tag structure
  void SetTimestamp (Time time);
  Time GetTimestamp (void) const;
  
  void Print (std::ostream &os) const;

private:
  Time m_timestamp;

  // end class TimestampTag
};

//----------------------------------------------------------------------

TypeId 
TimestampTag::GetTypeId (void)
{
  static TypeId tid = TypeId ("TimestampTag")
    .SetParent<Tag> ()
    .AddConstructor<TimestampTag> ()
    .AddAttribute ("Timestamp",
                   "Some momentous point in time!",
                   EmptyAttributeValue (),
                   MakeTimeAccessor (&TimestampTag::GetTimestamp),
                   MakeTimeChecker ())
  ;
  return tid;
}
TypeId 
TimestampTag::GetInstanceTypeId (void) const
{
  return GetTypeId ();
}

uint32_t 
TimestampTag::GetSerializedSize (void) const
{
  return 8;
}
void 
TimestampTag::Serialize (TagBuffer i) const
{
  int64_t t = m_timestamp.GetNanoSeconds ();
  i.Write ((const uint8_t *)&t, 8);
}
void 
TimestampTag::Deserialize (TagBuffer i)
{
  int64_t t;
  i.Read ((uint8_t *)&t, 8);
  m_timestamp = NanoSeconds (t);
}

void
TimestampTag::SetTimestamp (Time time)
{
  m_timestamp = time;
}
Time
TimestampTag::GetTimestamp (void) const
{
  return m_timestamp;
}

void 
TimestampTag::Print (std::ostream &os) const
{
  os << "t=" << m_timestamp;
}


// Special class of application that chnages the data rate every second according to a log-normal distribution

class MyApp : public Application 
{
public:

  MyApp ();
  virtual ~MyApp();

  void Setup (Ptr<Socket> socket, Address address, uint32_t packetSize, double mu, double sigma, std::vector<double> rate_mean, double rate_max, int Number_of_Intervals, double startTime, double endTime);
  
private:
  virtual void StartApplication (void);
  virtual void StopApplication (void);

  void ScheduleTx (void);
  void ScheduleStop (void);
  void SendPacket (void);
  void ChangeDataRate(void);

  Ptr<Socket>     m_socket;
  Address         m_peer;
  uint32_t        m_packetSize;
  DataRate        m_dataRate;
  double          m_mu;
  double          m_sigma;
  std::vector<double>          m_rate_mean;
  double 	  m_rate_max;
  int 		  m_Number_of_Intervals;
  EventId         m_sendEvent;
  EventId         m_changeRate;
  bool            m_running;
  uint32_t        m_packetsSent;
  double             m_startTime;
  double             m_endTime;
};

MyApp::MyApp ()
  : m_socket (0), 
    m_peer (), 
    m_packetSize (0),  
    m_dataRate (0), 
    m_mu (0),
    m_sigma (0),
    m_rate_mean (0),
    m_rate_max (0),
    m_Number_of_Intervals (0),
    m_sendEvent (), 
    m_changeRate (), 
    m_running (false), 
    m_packetsSent (0),
    m_startTime (0),
    m_endTime (0)
{
}

MyApp::~MyApp()
{
  m_socket = 0;
}

void
MyApp::Setup (Ptr<Socket> socket, Address address, uint32_t packetSize, double mu, double sigma, std::vector<double> rate_mean, double rate_max, int Number_of_Intervals, double startTime, double endTime)
{
  m_socket = socket;
  m_peer = address;
  m_packetSize = packetSize;
  m_dataRate = 0;
  m_mu = mu;
  m_sigma = sigma;
  m_rate_mean = rate_mean;
  m_Number_of_Intervals = Number_of_Intervals;
  m_startTime = startTime;
  m_endTime = endTime;
  m_rate_max = rate_max;
}

void
MyApp::StartApplication (void)
{
  nbAppRunning += 1;
  m_running = true;
  m_packetsSent = 0;
  m_socket->Bind ();
  m_socket->Connect (m_peer);
  ChangeDataRate();
  if(Pattern != "random"){
  	for (int s = 1; s < m_Number_of_Intervals; s++){ // remove these lines with random
  		m_changeRate = Simulator::Schedule (Seconds (SimStart-m_startTime + s * Interval_Adjustment), &MyApp::ChangeDataRate, this);
  	}
  }
  SendPacket ();
}

void 
MyApp::StopApplication (void)
{
  nbAppRunning -= 1;
  m_running = false;
  //std::cout << "Application closed, number of packets transmitted = " << m_packetsSent << " , still running = " << nbAppRunning << std::endl;

  if (m_sendEvent.IsRunning ())
    {
      Simulator::Cancel (m_sendEvent);
    }
    
  if (m_changeRate.IsRunning ())
    {
      Simulator::Cancel (m_changeRate);
    }

  if (m_socket)
    {
      m_socket->Close ();
    }
}

void 
MyApp::SendPacket (void)
{
  Ptr<Packet> packet = Create<Packet> (m_packetSize);
  TimestampTag timestamp;
  timestamp.SetTimestamp (Simulator::Now ());
  packet->AddByteTag (timestamp);
  
  m_socket->Send (packet);
  //std::cout << "Packet sent of size: " << packet-> GetSize()<< ", Available Space = " << m_socket -> GetTxAvailable() << std::endl;
  ScheduleTx ();

}

void 
MyApp::ScheduleTx (void)
{
  //std::cout << "Scheduling Tx " << m_running << " " << m_dataRate.GetBitRate ()<< std::endl;
  if (m_running)
    {
      Time tNext (Seconds (m_packetSize * 8 / static_cast<double> (m_dataRate.GetBitRate ())));
      m_sendEvent = Simulator::Schedule (tNext, &MyApp::SendPacket, this);
    }
}

void 
MyApp::ScheduleStop (void)
{
  if (m_running)
    {
      //Time tNext (Seconds (m_packetSize * 8 / static_cast<double> (m_dataRate.GetBitRate ())));
      Time tNext (Seconds (2.0));
      m_sendEvent = Simulator::Schedule (tNext, &MyApp::StopApplication, this);
    }
}

void 
MyApp::ChangeDataRate (void){
  //Ptr<LogNormalRandomVariable> x = CreateObject<LogNormalRandomVariable> ();
  //x->SetAttribute ("Mu", DoubleValue (m_mu));
  //x->SetAttribute ("Sigma", DoubleValue (m_sigma));
  
  Time timeNow = Simulator::Now();
  int index = int((timeNow.GetDouble()/1000000000 - SimStart)/Interval_Adjustment);
  //std::cout << index << std::endl;
  if (index < 0){index = 0;}
  double rate = m_rate_max * m_rate_mean[index];
  m_dataRate = DataRate (std::to_string(rate) + "Mbps");
  //std::cout << m_startTime - SimStart << " " << (m_startTime - SimStart)/(Interval_Adjustment) << " " << (int((m_startTime - SimStart)/(Interval_Adjustment))) << std::endl;
  //std::cout << "Data rate = " << rate << " Mbps at: "<< Simulator::Now() << std::endl;
}

void packetsRx (int i, Ptr<const Packet> packet,  Ptr< Ipv4 > ipv4, uint32_t interface)
{
  Ipv4Header ipv4h;
  //TcpHeader tcph;
  packet->PeekHeader (ipv4h);
  //packet->PeekHeader (tcph);
  
  //std::cout << "Packet of size: "<< packet -> GetSize() <<" received at: " << Simulator::Now() << ", at node: " << i << " "<< tcph << std::endl;
  	
  
  // Update node load
  load[i] += packet -> GetSize();
  NodeLoad[i] += packet -> GetSize();
  HeatR[RegionIds[i]-1] += packet -> GetSize() * 8/ Gscale; // in Gb
  
  // Update link load
  Ptr< NetDevice > ndI = ipv4 -> GetNetDevice (interface);
  for (int j = 0; j < totLinks; j++){
  	NetDeviceContainer ndT = ndc[j];
  	if (ndI == ndT.Get(0)){
		linkLoad[j] += packet -> GetSize();
  		break;
  	}else if (ndI == ndT.Get(1)){
  		linkLoad[j + totLinks] += packet -> GetSize();
  		break;
  	}
  } 
  
  // check if packet reached its destination
  int nInterfaces = nodes.Get(i)->GetNDevices(); 
  bool match = false;
  for (int index=0; index< nInterfaces; index++){
  Ipv4InterfaceAddress address_check = ipv4 -> GetAddress (index, 0);
  	if (ipv4h.GetDestination() == address_check.GetLocal()){
 	match = true;
  	break;
  	}
  }


  if (match){
  TimestampTag timestamp;
  if (packet->FindFirstMatchingByteTag (timestamp)) {
        Time tx = timestamp.GetTimestamp ();
        Time delay = Simulator::Now () - tx;
        double hopCount = 65 - (unsigned)ipv4h.GetTtl ();
        //std::cout << "At node: " << i << ", Delay = " << delay << ", Hop count = " << hopCount<< ", at time: " << Simulator::Now() << " of size: " << packet -> GetSize() <<std::endl;
        hopCount_array.push_back(hopCount);
        delay_array.push_back(delay.GetDouble());
        }
  }
}
  
void packetsTx (int i, Ptr<const Packet> packet,  Ptr< Ipv4 > ipv4, uint32_t interface)
{
  Ipv4Header ipv4h;
  //TcpHeader tcph;
  packet->PeekHeader (ipv4h);
  //packet->PeekHeader (tcph);
  
  //std::cout << "Packet of size : " << packet -> GetSize() <<" transmitted at: " << Simulator::Now() << ", from node: " << i << " " << tcph << std::endl;

  // check if this is the sending node
  int nInterfaces = nodes.Get(i)->GetNDevices(); 
  bool match = false;
  for (int index=0; index< nInterfaces; index++){
  Ipv4InterfaceAddress address_check = ipv4 -> GetAddress (index, 0);
  	if (ipv4h.GetSource() == address_check.GetLocal()){
 	match = true;
  	break;
  	}
  }
  
  if (match){
  	//std::cout << "This is a new packet!" << std::endl;
  	//std::cout << "Packet of size : " << packet -> GetSize() <<" generated at: " << Simulator::Now() << ", from node: " << i << " " << tcph << std::endl;
  	// Update node load:
  	load[i] += packet -> GetSize();
  	NodeLoad[i] += packet -> GetSize();
  	HeatR[RegionIds[i]-1] += packet -> GetSize() * 8 / Gscale; //in Gb
  	Throughput += packet -> GetSize(); 
  	if (packet -> GetSize() == 65535){
  		ThroughputD += packet -> GetSize();
  	}else if (packet -> GetSize() == 6000){
  		ThroughputR += packet -> GetSize();
  	}

  	for (uint32_t index = 0; index < Possible_size.size(); index++){
  		if(int(packet -> GetSize()) == Possible_size[index]){
  			Nb_Gen_perSize[index] +=1;
  		}  	
  	}
  	
  	//std::cout << ipv4h << std::endl;
  	//std::cout << tcph << std::endl;
  	/*
  	if (packet -> GetSize() != 56 and packet -> GetSize() != 52 and packet -> GetSize() != 65535 and packet -> GetSize() != 1500)
  	{
  		std::cout << "Packet of UNDESIRABLE size generated: " << packet -> GetSize() << std::endl;
  	} 
  	*/ 	
  }
  
}

void packetRxD(Ptr<const Packet> packet){
  //l0 += 1;
  std::cout << "Packet " << packet -> GetUid() << " of size: "<< packet -> GetSize() <<" received at: " << Simulator::Now() << std::endl;

}
void packetTxD(Ptr<const Packet> packet){
  //l0 -= 1;
  std::cout << "Packet " << packet -> GetUid() << " of size: "<< packet -> GetSize() <<" transmitted at: " << Simulator::Now()  << std::endl;

}

void trackChannel(Ptr< const Packet > packet, Ptr< NetDevice > txDevice, Ptr< NetDevice > rxDevice, Time duration, Time lastBitTime){
  std::cout << "Packet " << packet -> GetUid() << " of size: "<< packet -> GetSize() << ", duration = " << duration.GetDouble() << std::endl;
}


void packetsDrop (const Ipv4Header &header, Ptr< const Packet > packet,  Ipv4L3Protocol::DropReason reason, Ptr< Ipv4 > ipv4, uint32_t interface)
{
  dropped_packets += 1;
  std::cout<< "Packet dropped! Reason: " << reason << " , TTL = " << (unsigned)header.GetTtl () << std::endl;
}

void packetsDrop2 (Ptr< const Packet > packet)
{
  std::cout<< "Packet dropped at physical layer!" << std::endl;
}

void AdjustStart(){
  hopCount_array.clear();
  delay_array.clear();
  Throughput = 0;
  
  for(int i = 0; i < totNodes; i++){
  	load[i] = 0;
  	load_prev[i] = 0;
  	load_per_int[i] = 0;
  	NodeLoad[i] = 0;
  	Embedded_EC_perNode[i] = 0;
  	Operational_EC_perNode[i] = 0;
  	Embedded_CE_perNode[i] = 0;
  	Operational_CE_perNode[i] = 0;
  }
  
  for(int i = 0; i < 2 * totLinks; i++){
  	linkLoad[i] = 0;
  }
  
  for(int i = 0; i < totRegions; i++){
  	HeatR[i] = 0;
  } 
}


void UpdateThpt_MxLink(){
  Throughput = Throughput * 8 / 1000000000 / Interval_Adjustment; //Result is in Gbps
  Throughput_per_Interval.push_back(Throughput);
  //std::cout << "Throughput in this interval = " << Throughput << " Gbps" << std::endl;
  Throughput = 0;
  
  double max = 0;
  double total = 0;
  double Cap = 0;
  for (int i = 0; i < 2 * totLinks; i++){
  	double w;
  	if (i < totLinks){
  		w = weight[i];
  	}else{
  		w = weight[i - totLinks];
  	}
  	double val = linkLoad[i] * 8 / Interval_Adjustment / w / Gscale;  // /2
  	total += linkLoad[i] * 8 / Interval_Adjustment;
  	Cap += w * Gscale;  // *2
  	
  	if (val> max){
  		max = val;
  	}
  	linkLoad[i] = 0;
  }
  Max_linkLoad.push_back(max * 100);

}

void UpdateEnergy_Carbon(std::vector<double> CurrentForecast, int interval){

   for (int i=0; i < totNodes; i++){ 
      load_per_int[i] = load[i] - load_prev[i];
      // Calculate the new power consumption and carbon emissions per interval:
      double EC = (lamda[NodeLevel[i]-1] * load_per_int[i] * 8 / 1000000) / Interval_Adjustment ; // energy in Ws
      double CE = CurrentForecast[RegionIds[i]-1] * EC /3.6; // Carbon intensity is in ugCO2/Kwh
      // Update Global Variables
      EnergyConsumption += EC; // power in w
      CarbonEmissions += CE;
      // Update the per interval variables
      EnergyConsumption_per_intervalD[interval] += EC;
      CarbonEmissions_per_intervalD[interval] += CE;
      double EC2 = FixedNPower[NodeLevel[i]-1]; // in Ws /s
      double CE2 = CurrentForecast[RegionIds[i]-1] * EC2 / 3.6; // in ugCO2/Kwh /s
      EnergyConsumption_per_intervalF[interval] += EC2;
      CarbonEmissions_per_intervalF[interval] += CE2;
      // Update the per Node variables
      Embedded_EC_perNode[i] += EC2;
      Operational_EC_perNode[i] += EC;
      Embedded_CE_perNode[i] += CE2;
      Operational_CE_perNode[i] += CE;
   }
   load_prev = load;
   
   for (int i=0; i < totLinks; i++){ 
      	double EC = PortPower * 2; // power in Ws /s
      	double CE = (CurrentForecast[LinkRegions[i][0]-1] + CurrentForecast[LinkRegions[i][1]-1]) * PortPower / 3.6; 
      	EnergyConsumption_per_intervalP[interval] += EC;
        CarbonEmissions_per_intervalP[interval] += CE;
      	// Update the total energy consumption and carbon emissions
      	EnergyConsumption += EC; // power in w
        CarbonEmissions += CE;
   }
   
}

// Special Function to change the path costs:
void UpdatePathCosts(std::vector<double> CurrentForecast, int Scenario){
       std::cout << "Changing path costs! "<< std::endl;
       uint32_t nInterfaces = 0;
       for (int i=0; i < totNodes; i++){
       Ptr<Node> node = nodes.Get(i); 
       nInterfaces = nodes.Get(i)->GetNDevices();
       double cost = 1; // The cost is of type uint32_t which is ranged from 0 to 65535 --> should make sure of the range
       if (Scenario == 2){ // Incremental dynamic power per unit traffic
   	   cost = DPRating[NodeLevel[i]-1];
       } else if (Scenario == 3){ // Carbon Intensity
   	   cost = 1 + CurrentForecast[RegionIds[i]-1];
       } else if (Scenario == 4){ // Carbon Intensity + Incremental dynamic power per unit traffic
   	   cost = 1 + (CurrentForecast[RegionIds[i]-1]) * DPRating[NodeLevel[i]-1] / 10;
       } else if (Scenario == 5){ // Typical Power + Carbon Intensity
       	   cost = 1 + (CurrentForecast[RegionIds[i]-1]) * TypPower[NodeLevel[i]-1] / 10;
       } else if (Scenario == 6){ // Energy Labelling + Carbon Intensity
       	   cost = 1 + (CurrentForecast[RegionIds[i]-1]) * ERating[NodeLevel[i]-1];
       } else if (Scenario == 7){ // Typical Power
       	   cost = TypPower[NodeLevel[i]-1] / 10;
       } else if (Scenario == 8){ // Energy Labelling
       	   cost = ERating[NodeLevel[i]-1];
       } else if (Scenario == 9){ // Carbon Emissions
       	   cost = 1 + (CurrentForecast[RegionIds[i]-1]) * (FixedNPower[NodeLevel[i]-1] + (load[i]-load_prev[i]) * 8 /1000000 * lamda[NodeLevel[i]-1]) /10 ;
       	   //cost = 1 + (CurrentForecast[RegionIds[i]-1]) * (FixedNPower[NodeLevel[i]-1]) /10 ;
       	   //std::cout << cost << std::endl;
       } 
       
       
       for (uint32_t j=1; j < nInterfaces; j++){
            node -> GetObject<Ipv4> () -> SetMetric (j, cost);
            } 
       }
       
       // Flip the cost metric for every pair of interfaces
       // A router sends the packet through the interface that has the lower cost BUT that cost is the one set by the neighbor node -> that is why we flip the cost of interfaces.
       for (int i=0; i < totLinks; i++){
       uint16_t cost0 = ipic[i].Get (0).first -> GetObject<Ipv4L3Protocol> ()->GetInterface (ipic[i].Get (0).second) -> GetMetric();
      uint16_t cost1 = ipic[i].Get (1).first -> GetObject<Ipv4L3Protocol> ()->GetInterface (ipic[i].Get (1).second) -> GetMetric();
      ipic[i].SetMetric (0, cost1);
      ipic[i].SetMetric (1, cost0);
       }  
}

// Special Function to get the forecast:
std::vector<std::vector<double>> GetForecast (std::string Season)
{
  std::string FileName (dir + Topology +"Files/" + Season + ".txt");
  std::ifstream topgen;
  topgen.open (FileName);
  std::vector<std::vector<double>> Forecast;

  if ( !topgen.is_open () )
    {
      NS_LOG_WARN ("Forecast file object is not open, check file name and permissions");
      return Forecast;
    }

  std::istringstream lineBuffer;
  std::string line;
  
  int value;

  if (Pattern == "random"){
	  for (int i = 0; i < Number_of_Intervals && !topgen.eof (); i++){
	  	getline (topgen,line);
	      	lineBuffer.clear ();
	      	lineBuffer.str (line);
	      	std::vector <double> f;
	      	for(int j = 0; j < totRegions; j++){
	      		lineBuffer >> value;
	      		f.push_back(value);
	  	}
	   	Forecast.push_back(f); 
	   }
  }else{
  	for (int i = 0; i < 16 + Number_of_Intervals && !topgen.eof (); i++){
	  	getline (topgen,line);
	  	if (i >= 16){
		      	lineBuffer.clear ();
		      	lineBuffer.str (line);
		      	std::vector <double> f;
		      	for(int j = 0; j < totRegions; j++){
		      		lineBuffer >> value;
		      		f.push_back(value);
		  	}
		   	Forecast.push_back(f); 
	   		}
  		}
  	}	  
   
   topgen.close ();

   return Forecast;
}


// Special Function to get the rate mean associated with the topology:
std::vector<double> GetRateMean ()
{
  std::string FileName (dir + Topology +"Files/RateMean.txt");
  std::ifstream topgen;
  topgen.open (FileName);
  std::vector<double> RateMean;

  if ( !topgen.is_open () )
    {
      NS_LOG_WARN ("Rate Mean file object is not open, check file name and permissions");
      return RateMean;
    }

  std::istringstream lineBuffer;
  std::string line;
  
  double value;
  getline (topgen,line);
  lineBuffer.clear ();
  lineBuffer.str (line);
  if (Pattern == "random"){
	  for(int j = 0; j < Number_of_Intervals; j++){
	      	lineBuffer >> value;
	      	RateMean.push_back(value);
	  }
  }else{
  	for(int j = 0; j < 16 + Number_of_Intervals; j++){
	      	lineBuffer >> value;
	      	if (j >= 16){
	      		RateMean.push_back(value);
	      	}
	  }
  }
  
  topgen.close ();
  return RateMean;
}

// Special Function to read the topolgy based on Inet Format
void CreateTopology (std::string FileName)
{
  std::ifstream topgen;
  topgen.open (FileName);

  if ( !topgen.is_open () )
    {
      NS_LOG_WARN ("Inet topology file object is not open, check file name and permissions");
    }

  std::istringstream lineBuffer;
  std::string line;

  getline (topgen,line);
  lineBuffer.str (line);
  lineBuffer >> totNodes;
  lineBuffer >> totLinks;
  lineBuffer >> totRegions;
  
  nodes.Create(totNodes);
  std::vector <int> NodeLinks;

  for (int i = 0; i < totNodes && !topgen.eof (); i++)
    {
      getline (topgen,line);
      lineBuffer.clear ();
      lineBuffer.str (line);
      
      std::string node_id;
      int NodeRegionId;
      int node_level;
      
      lineBuffer >> node_id;
      lineBuffer >> NodeRegionId;
      lineBuffer >> node_level;
      RegionIds.push_back(NodeRegionId);
      NodeLevel.push_back(node_level);
      AgRate.push_back(0);
      AllNodes.push_back(i);
      NodeLinks.push_back(0);
      if (node_level == 1){
      	CoreNodes.push_back(i);
      }
      if (node_level == 2){
      	MetroNodes.push_back(i); 
      }
      if (node_level == 3){
      	Tier1Nodes.push_back(i);
      }
      
      //std::cout << "Node: " << node_id << ", Region: " << NodeRegionId << ", level: " << node_level << std::endl;
    }

  for (int i = 0; i < totLinks && !topgen.eof (); i++)
    {
      getline (topgen,line);
      lineBuffer.clear ();
      lineBuffer.str (line);
      
      int fromNode;
      int toNode;
      double r;
      double d;
      
      lineBuffer >> fromNode;
      lineBuffer >> toNode;
      lineBuffer >> r;
      lineBuffer >> d;
      
      nc.push_back(NodeContainer (nodes.Get(fromNode), nodes.Get(toNode)));
      weight.push_back(r /Rscale); // in Gbps
      delay_perLink.push_back(std::to_string(5 * d)); // 5 ns per m // default time unit is ns
      AgRate[fromNode] += r / Rscale;
      AgRate[toNode] += r / Rscale;
      LinkRegions.push_back({RegionIds[fromNode], RegionIds[toNode]});
      NodeLinks[fromNode] += 1;
      NodeLinks[toNode] += 1;
      //std::cout << "from Node: " << fromNode << " (level: " << NodeLevel[fromNode] << "), To Node: " << toNode << " (level: " << NodeLevel[toNode] << "), Rate: " << r << " Gbps, Delay: " << 5 * d << " ns" << std::endl;
    }
  
  topgen.close ();

}


// Special Function to read the topolgy parameters
void ReadParameters (std::string FileName)
{
  std::ifstream topgen;
  topgen.open (FileName);

  if ( !topgen.is_open () )
    {
      NS_LOG_WARN ("Parameters file object is not open, check file name and permissions");
    }

  std::istringstream lineBuffer;
  std::string line;
  
  getline (topgen,line);
  lineBuffer.clear ();
  lineBuffer.str (line);
  std::string id;
  lineBuffer >> id;
  lineBuffer >> id;
  lineBuffer >> Pattern;
  
  getline (topgen,line);
  lineBuffer.clear ();
  lineBuffer.str (line);
  int NodeHierarchyLevels;
  lineBuffer >> id;
  lineBuffer >> id;
  lineBuffer >> NodeHierarchyLevels;
  
  std::vector<float> list = {};
  
  for (int i = 0; i < 13 && !topgen.eof (); i++)
    {
      getline (topgen,line);
      lineBuffer.clear ();
      lineBuffer.str (line);
      std::string id;
      float val;
      lineBuffer >> id;
      lineBuffer >> id;
      lineBuffer >> val;
      list.push_back(val);
      if (i >= 8 and NodeHierarchyLevels > 1){
      	for (int j = 1; j < NodeHierarchyLevels; j++){
      		lineBuffer >> val;
      		list.push_back(val);
      	}
      }
    }
    
    int i = 0;
    Sec_per_Interval = list[i++]; 
    Interval_Adjustment = Sec_per_Interval;
    Number_of_Intervals = int(list[i++]); 
    rate_max = list[i++]; 
    nb_Applications = int(list[i++]); 
    SimStart = list[i++]; 
    CustomPacketSize = int(list[i++]); 
    offset = list[i++]; 
    PortPower = list[i++]; 
    FixedNPower = {list[i++]}; 
    if (NodeHierarchyLevels > 1){
      	for (int j = 1; j < NodeHierarchyLevels; j++){
      		FixedNPower.push_back(list[i++]); 
      	}
      }
    ERating = {list[i++]}; 
    if (NodeHierarchyLevels > 1){
      	for (int j = 1; j < NodeHierarchyLevels; j++){
      		ERating.push_back(list[i++]); 
      	}
      }
    TypPower = {list[i++]}; 
    if (NodeHierarchyLevels > 1){
      	for (int j = 1; j < NodeHierarchyLevels; j++){
      		TypPower.push_back(list[i++]); 
      	}
      }
    lamda = {list[i++]}; 
    if (NodeHierarchyLevels > 1){
      	for (int j = 1; j < NodeHierarchyLevels; j++){
      		lamda.push_back(list[i++]); 
      	}
      }
    DPRating = {list[i++]}; 
    if (NodeHierarchyLevels > 1){
      	for (int j = 1; j < NodeHierarchyLevels; j++){
      		DPRating.push_back(list[i++]); 
      	}
      }
  topgen.close ();

}

// Special Function to build the traffix matrix in case of downstreaming
void CreateDownStreamPattern (std::string FileName)
{
  std::ifstream topgen;
  topgen.open (FileName);

  if ( !topgen.is_open () )
    {
      NS_LOG_WARN ("Inet topology file object is not open, check file name and permissions");
    }

  std::istringstream lineBuffer;
  std::string line;
  getline (topgen,line);

  for (int i = 0; i < totNodes && !topgen.eof (); i++)
    {
      getline (topgen,line);
      Neighbors.push_back({});
    }

  for (int i = 0; i < totLinks && !topgen.eof (); i++)
    {
      getline (topgen,line);
      lineBuffer.clear ();
      lineBuffer.str (line);
      
      int fromNode;
      int toNode;
      double r;
      double d;
      
      lineBuffer >> fromNode;
      lineBuffer >> toNode;
      lineBuffer >> r;
      lineBuffer >> d;
      Neighbors[fromNode].push_back(toNode);
      Neighbors[toNode].push_back(fromNode);
      
    }
    
    for (uint32_t i = 0; i < CoreNodes.size(); i++){
    	int node = CoreNodes[i];
    	int c = 0; //counter
    	for (uint32_t j = 0; j < Neighbors[node].size(); j++){
    		int n = Neighbors[node][j];
    		if (NodeLevel[n] == 2){
    			Flows.push_back({n, node});
    			c ++;
    			if (c == 2){break;}
    		}
    	}
    }
    for (uint32_t i = 0; i < MetroNodes.size(); i++){
    	int node = MetroNodes[i];
    	int c = 0; //counter
    	for (uint32_t j = 0; j < Neighbors[node].size(); j++){
    		int n = Neighbors[node][j];
    		if (NodeLevel[n] == 2){
    			Flows.push_back({n, node});
    			c ++;
    			if (c == 1){break;}
    		}
    	}
    	if(c < 1){
    		for (uint32_t j = 0; j < Neighbors[node].size(); j++){
    			if (c == 1){break;}
    			int n = Neighbors[node][j];
    			std::vector<int> FurtherNeighbors = Neighbors[n];
    			for (uint32_t z = 0; z < FurtherNeighbors.size(); z++){
    				if (FurtherNeighbors[z] == node){continue;}
    				if (NodeLevel[FurtherNeighbors[z]] == 2){
    					Flows.push_back({FurtherNeighbors[z], node});
    					c ++;
    					if (c == 1){break;}
    				} 			
    			}
    			
    		}
    	}
    }
    
    for (uint32_t i = 0; i < Tier1Nodes.size(); i++){
    	int node = Tier1Nodes[i];
    	int c = 0; //counter
    	for (uint32_t j = 0; j < Neighbors[node].size(); j++){
    		int n = Neighbors[node][j];
    		if (NodeLevel[n] == 2){
    			Flows.push_back({n, node});
    			c ++;
    			if (c == 2){break;}
    		}
    	}
    	if(c < 2){
    		for (uint32_t j = 0; j < Neighbors[node].size(); j++){
    			if (c == 2){break;}
    			int n = Neighbors[node][j];
    			std::vector<int> FurtherNeighbors = Neighbors[n];
    			for (uint32_t z = 0; z < FurtherNeighbors.size(); z++){
    				if (FurtherNeighbors[z] == node){continue;}
    				if (NodeLevel[FurtherNeighbors[z]] == 2){
    					Flows.push_back({FurtherNeighbors[z], node});
    					c ++;
    					if (c == 2){break;}
    				} 			
    			}
    			
    		}
    	}
    }

  
  topgen.close ();

}

// ----------------------------------------------------------------------
// -- main
// ----------------------------------------------
int main (int argc, char *argv[])
{
  Config::SetDefault("ns3::TcpL4Protocol::SocketType", StringValue("ns3::TcpCubic"));
  std::string format ("Inet");
  
  // Set up command line parameters used to control the experiment.
  int Scenario = 1; //1 to 9 as explained in the changePathCost function
  std::string Season = "fall";
  CommandLine cmd (__FILE__);
  cmd.AddValue ("Topology", "BT or GEANT or ...", Topology);
  cmd.AddValue ("Scenario", "from 1 to 9", Scenario);
  cmd.AddValue ("Season", "fall, winter, spring, summer", Season);
  cmd.AddValue ("Directory", "/home/username/", dir);
  cmd.Parse (argc, argv);
  
  dir = dir + "ns-allinone-3.36.1/ns-3.36.1/scratch/CATE/";
  std::cout << "Topology = " << Topology << std::endl;
  std::cout << "Scenario = " << Scenario << std::endl;
  std::cout << "Season = " << Season << std::endl;
  
  std::string input;
  input = dir + Topology + "Files/parameters.txt";
  ReadParameters(input); // Adjust the simulation parameters based on the file located under the TOPOLOGYFiles directory.
  
  std::cout << "Pattern = " << Pattern << std::endl;
  
  double SimTimeSec = SimStart + Sec_per_Interval * Number_of_Intervals; // in Seconds.
  std::cout << "Simulation Time = " << SimTimeSec << std::endl;
  double ExtraTimeSec = 0; // can add for sink app and the overall simulation but carefull with indices
  
  input = dir + Topology + "Files/Topology.txt";
  
  
  //LogComponentEnableAll(LOG_LEVEL_INFO);
  Config::SetDefault ("ns3::Ipv4GlobalRouting::RespondToInterfaceEvents", BooleanValue (true));
  Config::SetDefault ("ns3::Ipv4GlobalRouting::RandomEcmpRouting", BooleanValue (true));
  
  // ------------------------------------------------------------
  // -- Read topology data.
  // --------------------------------------------

  // Pick a topology reader based in the requested format.
  CreateTopology (input);
 
  //Create network stacks, node containers and corresponding net devices and interfaces
  NS_LOG_INFO ("creating internet stack");
  
  InternetStackHelper stack;
  stack.Install (nodes);
  

  NS_LOG_INFO ("creating IPv4 addresses");
  Ipv4AddressHelper address;
  address.SetBase ("10.0.0.0", "255.255.255.252");


  NS_LOG_INFO ("creating net device containers");
  PointToPointHelper p2p;
  for (int i = 0; i < totLinks; i++)
    {
      p2p.SetChannelAttribute ("Delay", StringValue (delay_perLink[i] +"ns"));
      p2p.SetDeviceAttribute ("DataRate", StringValue (std::to_string(weight[i])+"Gbps")); 
      NetDeviceContainer ndcT = p2p.Install (nc[i]);
      ndcT.Get(0) -> SetMtu(CustomPacketSize); // Increase the packet size to increase the utilization of the system with faster simulation
      ndcT.Get(1) -> SetMtu(CustomPacketSize);
      ndc.push_back(ndcT);
    }
     
  

  NS_LOG_INFO ("creating IPv4 interfaces");
  for (int i = 0; i < totLinks; i++)
    {
      ipic.push_back(address.Assign (ndc[i]));
      address.NewNetwork ();
    }

  // Get the total Carbon Emissions forecast (interval: 30 min)
  std::vector<std::vector<double>> Forecast = GetForecast(Season);

  // Track the load at every node per second:
  for(int i = 0; i < totNodes; i++){
  	load.push_back(0);
  	load_prev.push_back(0);
  	load_per_int.push_back(0);
  	NodeLoad.push_back(0);
  	Embedded_EC_perNode.push_back(0);
  	Operational_EC_perNode.push_back(0);
  	Embedded_CE_perNode.push_back(0);
  	Operational_CE_perNode.push_back(0);
  }
  
  for(int i = 0; i < 2 * totLinks; i++){
  	linkLoad.push_back(0);
  }
  
  for(int i = 0; i < totRegions; i++){
  	HeatR.push_back(0);
  }
  
  // Track the energy consumption and carbon emissions per interval:
  for(int i = 0; i < Number_of_Intervals; i++){
  	EnergyConsumption_per_intervalD.push_back(0);
  	CarbonEmissions_per_intervalD.push_back(0);
  	EnergyConsumption_per_intervalP.push_back(0);
  	CarbonEmissions_per_intervalP.push_back(0);
  	EnergyConsumption_per_intervalF.push_back(0);
  	CarbonEmissions_per_intervalF.push_back(0);
  }


  // Update the path costs before computing the routing tables: (except for scenarios with utilization)
  if (Scenario != 1){
  	UpdatePathCosts(Forecast[0], Scenario);
  }
  
  
  // Compute the initial routing tables
  std::cout << "Before Computing Routing Tables \n";
  Ipv4GlobalRoutingHelper::PopulateRoutingTables (); 
  std::cout << "After Computing Routing Tables \n";
  
  
  // Generate Log-normal traffic at the lowest level nodes
  NS_LOG_INFO ("Create Applications.");

  //Uniform random variable to get a random server node
  //Ptr<UniformRandomVariable> unifRandom = CreateObject<UniformRandomVariable> ();
  //unifRandom->SetAttribute ("Min", DoubleValue (0));
  //unifRandom->SetAttribute ("Max", DoubleValue (totNodes - 1));
  

  double mu = 0;
  double sigma = 0.25;
  std::vector<double> rate_mean = GetRateMean();

  
  uint32_t payloadSize = CustomPacketSize - 52; // in bytes // this is the maximum allowed payloadsize to comply with MTU
  Config::SetDefault ("ns3::TcpSocket::SegmentSize", UintegerValue (payloadSize)); 
  
  ApplicationContainer apps;
  int randomServerNumber = 0;
  
  // Creating sink at all nodes
  uint16_t port = 7;    
  Address localAddress (InetSocketAddress (Ipv4Address::GetAny (), port));
  PacketSinkHelper packetSinkHelper ("ns3::TcpSocketFactory", localAddress);
  ApplicationContainer sinkApp;
  
  for (int i = 0; i < totNodes; i++){
  sinkApp.Add(packetSinkHelper.Install (nodes.Get (i)));
  }
  
  sinkApp.Start (Seconds (0));
  sinkApp.Stop (Seconds (SimTimeSec + ExtraTimeSec));
  
  
  // Choose which nodes are among the traffic generating nodes
  std::vector<int> TrafficNodes = AllNodes;
  int totTrafficNodes = TrafficNodes.size();
  std::cout << "Traffic nodes size = "<< totTrafficNodes << std::endl;
  // Creating sockets at all nodes
  Ptr<Socket>* socket = new Ptr<Socket>[totTrafficNodes];
  
  if (Pattern == "random"){
  double spacing = 0.000005; //comparable to the transmission time of a packet
  std::srand(0);  // Initialize random number generator.
  for (int counter = 0; counter< int((SimTimeSec-offset) / spacing); counter++){ 
  	for (int a = 0; a < nb_Applications; a++){
	  int src = rand() % totNodes;
	  int dst = 0;
	  do{
	  	dst = rand() % totNodes;
	  }while(dst == src);
	  Ptr<Node> randomServerNode = nodes.Get (dst);
	  Ptr<Ipv4> ipv4Server = randomServerNode->GetObject<Ipv4> ();
	  Ipv4InterfaceAddress iaddrServer = ipv4Server->GetAddress (1,0);
	  Ipv4Address ipv4AddrServer = iaddrServer.GetLocal ();
	  InetSocketAddress rmt (ipv4AddrServer, port);
	  
	  socket[src] = Socket::CreateSocket (nodes.Get (src), TcpSocketFactory::GetTypeId ());
	  Ptr<MyApp> app = CreateObject<MyApp> ();
	  double startTime = offset + counter * spacing; 
	  double endTime = startTime + 0.055; // time to stop the application, can be changed
	  app->Setup (socket[src], rmt, payloadSize, mu, sigma, rate_mean, rate_max, Number_of_Intervals, startTime, endTime);
	  nodes.Get (src)->AddApplication (app);
	  app->SetStartTime (Seconds (startTime));
	  app->SetStopTime (Seconds(endTime));
  	}
  
  }
  }else if (Pattern == "downstreaming"){
  	double counter = 0;
  	CreateDownStreamPattern (input);
  	std::cout << "FLows number = " << Flows.size() << std::endl;
  	for (uint32_t i = 0; i< Flows.size(); i++){
  		randomServerNumber = Flows[i][1];
		Ptr<Node> randomServerNode = nodes.Get (randomServerNumber);
	  	Ptr<Ipv4> ipv4Server = randomServerNode->GetObject<Ipv4> ();
	  	Ipv4InterfaceAddress iaddrServer = ipv4Server->GetAddress (1,0);
	  	Ipv4Address ipv4AddrServer = iaddrServer.GetLocal ();
	  	InetSocketAddress rmt (ipv4AddrServer, port);
	  	for (int a = 0; a < nb_Applications; a++){
			socket[Flows[i][0]] = Socket::CreateSocket (nodes.Get (Flows[i][0]), TcpSocketFactory::GetTypeId ());
			Ptr<MyApp> app = CreateObject<MyApp> ();
			double starttime = offset + counter * 0.000005; //comparable to the transmission time of a packet
			counter += 1;
	  		app->Setup (socket[Flows[i][0]], rmt, payloadSize, mu, sigma, rate_mean, rate_max, Number_of_Intervals, starttime, SimTimeSec);	  	
	  		nodes.Get (Flows[i][0])->AddApplication (app);
	  		app->SetStartTime (Seconds (starttime));	
	  	}
  	}
  }else{ // mixed random and downstreaming
  
  	// Downstreaming Traffic part
  	rate_max = rate_max * 0.7; //you can scale up or down
  	
  	double counter = 0;
  	CreateDownStreamPattern (input);
  	std::cout << "FLows number = " << Flows.size() << std::endl;
  	for (uint32_t i = 0; i< Flows.size(); i++){
  		randomServerNumber = Flows[i][1];
		Ptr<Node> randomServerNode = nodes.Get (randomServerNumber);
	  	Ptr<Ipv4> ipv4Server = randomServerNode->GetObject<Ipv4> ();
	  	Ipv4InterfaceAddress iaddrServer = ipv4Server->GetAddress (1,0);
	  	Ipv4Address ipv4AddrServer = iaddrServer.GetLocal ();
	  	InetSocketAddress rmt (ipv4AddrServer, port);
	  	for (int a = 0; a < nb_Applications; a++){
			socket[Flows[i][0]] = Socket::CreateSocket (nodes.Get (Flows[i][0]), TcpSocketFactory::GetTypeId ());
			Ptr<MyApp> app = CreateObject<MyApp> ();
			double starttime = offset + counter * 0.000005; //comparable to the transmission time of a packet
			counter += 1;
	  		app->Setup (socket[Flows[i][0]], rmt, payloadSize, mu, sigma, rate_mean, rate_max, Number_of_Intervals, starttime, SimTimeSec);	  	
	  		nodes.Get (Flows[i][0])->AddApplication (app);
	  		app->SetStartTime (Seconds (starttime));	
	  	}
  	}
  	
  	// Random traffic part 
  	
	rate_max = rate_max * 0.4; 
	nb_Applications = 5;
	CustomPacketSize = 6000;
    	double offset2 = 0.95;
    	
    	double spacing = 0.000005; 
  	std::srand(0);  // Initialize random number generator.
  	for (int counter = 0; counter< int((SimTimeSec-offset2) / spacing); counter++){ 
	  	for (int a = 0; a < nb_Applications; a++){
		  int src = rand() % totNodes;
		  int dst = 0;
		  do{
		  	dst = rand() % totNodes;
		  }while(dst == src);
		  Ptr<Node> randomServerNode = nodes.Get (dst);
		  Ptr<Ipv4> ipv4Server = randomServerNode->GetObject<Ipv4> ();
		  Ipv4InterfaceAddress iaddrServer = ipv4Server->GetAddress (1,0);
		  Ipv4Address ipv4AddrServer = iaddrServer.GetLocal ();
		  InetSocketAddress rmt (ipv4AddrServer, port);
		  
		  socket[src] = Socket::CreateSocket (nodes.Get (src), TcpSocketFactory::GetTypeId ());
		  Ptr<MyApp> app = CreateObject<MyApp> ();
		  double startTime = offset2 + counter * spacing; //comparable to the transmission time of a packet
		  double endTime = startTime + 0.055; 
		  app->Setup (socket[src], rmt, CustomPacketSize - 52, mu, sigma, rate_mean, rate_max, Number_of_Intervals, startTime, endTime);
		  nodes.Get (src)->AddApplication (app);
		  app->SetStartTime (Seconds (startTime));
		  app->SetStopTime (Seconds(endTime));
	  	}
  	}
  }	
  
  
  // 
  for(int i = 0; i < totNodes; i++){
  std::string str1 = "/NodeList/" + std::to_string(i) + "/$ns3::Ipv4L3Protocol/Tx";
  std::string str2 = "/NodeList/" + std::to_string(i) + "/$ns3::Ipv4L3Protocol/Rx";
  Config::ConnectWithoutContext(str1, MakeBoundCallback(&packetsTx,i));
  Config::ConnectWithoutContext(str2, MakeBoundCallback(&packetsRx,i));
  }
  Config::ConnectWithoutContext("/NodeList/*/$ns3::Ipv4L3Protocol/Drop", MakeCallback(&packetsDrop));
  
    
  // Reset at 1 second which is the actual start time // the flows started at offset time so that at 1s the system is at a normal state (without this, the first interval of simulation will start with empty buffers and will take time to be at a normal state)
  if (Scenario == 9){
  	Simulator::Schedule (Seconds (SimStart), &UpdatePathCosts, Forecast[0], Scenario);
  	Simulator::Schedule (Seconds (SimStart), Ipv4GlobalRoutingHelper::RecomputeRoutingTables);
  }
  Simulator::Schedule (Seconds (SimStart), &AdjustStart);
  
  // Update the path costs for every interval only for scenarios with carbon intensity considered:
  if (Scenario == 3 or Scenario == 4 or Scenario == 5 or Scenario == 6 or Scenario == 9){
  	for (int interval = 1; interval < Number_of_Intervals; interval++){
     		Simulator::Schedule (Seconds (SimStart + Interval_Adjustment * interval), &UpdatePathCosts, Forecast[interval], Scenario);
     		Simulator::Schedule (Seconds (SimStart + Interval_Adjustment * interval), Ipv4GlobalRoutingHelper::RecomputeRoutingTables);
  	}
  }
  
  // Update the energy consumption, the carbon emissions, the throughput and the max link load for every interval:
  for (int interval = 1; interval <= Number_of_Intervals; interval++){
  	Simulator::Schedule (Seconds (SimStart + Interval_Adjustment * interval), &UpdateEnergy_Carbon, Forecast[interval-1], interval-1);
  	Simulator::Schedule (Seconds (SimStart + Interval_Adjustment * interval), &UpdateThpt_MxLink);
  }
  
  std::cout << "Simulation started" << std::endl;
  // Run the Simulation
  NS_LOG_INFO ("Run Simulation.");
  Simulator::Stop (Seconds (SimTimeSec + ExtraTimeSec));   // Define the end time of simulation, necessary for flowmonitor, not used here
  Simulator::Run ();
  
  double APD = 0;
  double STD = 0;
  double TailLatency = 0;
  double avgHopCount = 0;
  double STH = 0;
  double TailHopCount = 0;
  // Read the statistics collected: APD, p99 Tail Latency, Energy Consumption and Carbon Emissions
  if (delay_array.size() !=0){
  	std::cout << "Delay array size = " << delay_array.size() << std::endl;
  	APD = std::accumulate(delay_array.begin(), delay_array.end(), 0.0) / delay_array.size();
  	std::cout << "Average Packet delay = " << APD / 1000000 << " ms" << "\n";
  	for(uint32_t i = 0; i < delay_array.size(); i++){
  		STD += std::pow(delay_array[i] - APD, 2);}
  	STD = std::sqrt(STD / (delay_array.size() - 1));
  	std::cout << "Standard Deviation of Delay = " << STD / 1000000 << " ms" << std::endl;
  	std::sort(delay_array.begin(), delay_array.end());
  	TailLatency = delay_array[int(0.99*(delay_array.size()-1))];
  	std::cout << "p99 Tail Latency = " << TailLatency / 1000000 << " ms \n";
  }
  if (hopCount_array.size() !=0){
  	std::cout << "Hop Count array size = " << hopCount_array.size() << std::endl;
  	avgHopCount = std::accumulate(hopCount_array.begin(), hopCount_array.end(), 0.0) / hopCount_array.size();
  	std::cout << "Average Hop Count = " << avgHopCount << "\n";
  	for(uint32_t i = 0; i < hopCount_array.size(); i++){
  		STH += std::pow(hopCount_array[i] - avgHopCount, 2);}
  	STH = std::sqrt(STH / (hopCount_array.size() - 1));
  	std::cout << "Standard Deviation of Hop Count = " << STH << std::endl;
  	std::sort(hopCount_array.begin(), hopCount_array.end());
  	TailHopCount = hopCount_array[int(0.99*(hopCount_array.size()-1))];
  	std::cout << "p99 Hop Count Tail = " << TailHopCount << "\n";
  }

  

  double total_link_capacity = std::accumulate(weight.begin(), weight.end(), 0.0);
  double total_linkLoad = std::accumulate(linkLoad.begin(), linkLoad.end(), 0.0);
  std::cout << "Energy Consumption = " << EnergyConsumption /3600000 << " KWh \n";
  std::cout << "Carbon Emissions = " << CarbonEmissions/1000000 << " gCO2 \n";
  std::cout << "Number of dropped packets = " << dropped_packets << std::endl;
  std::cout << "Throughput for downstreaming = " << ThroughputD * 8 /1000000000 << " Gbps" << std::endl;
  std::cout << "Throughput for random = " << ThroughputR * 8 /1000000000 << " Gbps" << std::endl;
  for (uint32_t index = 0; index < Possible_size.size(); index++){
  	std::cout << "For size " << Possible_size[index] << " nb = " << Nb_Gen_perSize[index] << std::endl;	 	
  }
 
  for (int interval = 0; interval < Number_of_Intervals; interval++){
  	std::cout << interval << " Throughput = " << Throughput_per_Interval[interval] << " Gbps, max link load = " << Max_linkLoad[interval] << std::endl;
  }

  
  //std::cout << std::accumulate(CarbonEmissions_per_intervalF.begin(), CarbonEmissions_per_intervalF.end(), 0.0) /1000000 * 30 * 60 << " " << std::accumulate(CarbonEmissions_per_intervalP.begin(), CarbonEmissions_per_intervalP.end(), 0.0) /1000000 * 30 * 60 << " " << std::accumulate(CarbonEmissions_per_intervalD.begin(), CarbonEmissions_per_intervalD.end(), 0.0) /1000000 * 30 * 60 << std::endl;
  
  // Write into new txt files: delay_array, hopCount_array, NodeLoad/AgRate, linkLoad/weight, energy consumption per interval, carbon emissions per interval.
  // report the global values
  std::string file1 = dir + Topology +"Files/Results/DelayValues" + Pattern + Season + std::to_string(Scenario) + ".txt";
  std::string file2 = dir + Topology +"Files/Results/HopCountValues" + Pattern + Season + std::to_string(Scenario) + ".txt";
  std::string file3 = dir + Topology +"Files/Results/Results" + Pattern + Season + std::to_string(Scenario) + ".txt";
  FILE *fpt;
  
  const char* str = file1.c_str();
  fpt = fopen(str, "w"); 
  
  for(uint32_t i = 0; i < delay_array.size(); i++){
  	fprintf(fpt,"%lf \n", delay_array[i]/1000000); // ms
  }
  fclose(fpt);
  
  str = file2.c_str();
  fpt = fopen(str, "w"); 
  
  for(uint32_t i = 0; i < delay_array.size(); i++){
  	fprintf(fpt,"%lf \n", hopCount_array[i]);
  }
  fclose(fpt);

  

  
  str = file3.c_str();
  fpt = fopen(str, "w"); 
  
  fprintf(fpt, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", CarbonEmissions/1000000, EnergyConsumption/3600000, APD / 1000000, STD / 1000000, TailLatency / 1000000, avgHopCount, STH, TailHopCount);
  
  
  for(int i=0; i < Number_of_Intervals; i++){
  	fprintf(fpt,"%lf \t", EnergyConsumption_per_intervalD[i]/3600000); // in KWh per second
  }
  fprintf(fpt,"\n");

  for(int i=0; i < Number_of_Intervals; i++){
  	fprintf(fpt,"%lf \t", CarbonEmissions_per_intervalD[i]/1000000); // in gCO2 per second
  }
  fprintf(fpt,"\n");

  for(int i=0; i < Number_of_Intervals; i++){
  	fprintf(fpt,"%lf \t", Throughput_per_Interval[i]); // in Gbps
  }
  fprintf(fpt,"\n");

  for(int i=0; i < Number_of_Intervals; i++){
  	fprintf(fpt,"%lf \t", Max_linkLoad[i]); // in %
  }
  fprintf(fpt,"\n");

  for(int i = 0; i < totRegions; i++){
  	fprintf(fpt,"%lf\t", HeatR[i] /(SimTimeSec-SimStart)); // expressed in Gbps throughput overall avg in all intervals
  }
  fprintf(fpt,"\n");  
  
  for(int i=0; i < Number_of_Intervals; i++){
  	fprintf(fpt,"%lf \t", EnergyConsumption_per_intervalF[i]/3600000); // in KWh per second
  }
  fprintf(fpt,"\n");

  for(int i=0; i < Number_of_Intervals; i++){
  	fprintf(fpt,"%lf \t", CarbonEmissions_per_intervalF[i]/1000000); // in gCO2 per second
  }
  fprintf(fpt,"\n");
  
  for(int i=0; i < Number_of_Intervals; i++){
  	fprintf(fpt,"%lf \t", EnergyConsumption_per_intervalP[i]/3600000); // in KWh per second
  }
  fprintf(fpt,"\n");

  for(int i=0; i < Number_of_Intervals; i++){
  	fprintf(fpt,"%lf \t", CarbonEmissions_per_intervalP[i]/1000000); // in gCO2 per second
  }
  fclose(fpt);

  Simulator::Destroy ();
  
  delete[] socket;

  return 0;

  // end main
}
