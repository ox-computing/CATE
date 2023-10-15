from priodict import priorityDictionary
    
    
def Dijkstra(ECMP, G,start,end=None):

    D = {}  # dictionary of final distances
    list_P = []
    P = {}  # dictionary of predecessors
    list_P.append(P)
    Q = priorityDictionary()
    Q[start] = 0
    ECMP_tuples1 = []
    ECMP_tuples2 = []
    for v in Q:
        D[v] = Q[v]
        if v == end: break
        for w in G[v]:
            vwLength = D[v] + G[v][w]
            if w in D:
                if vwLength < D[w]:
                    raise ValueError ("Dijkstra: found better path to already-final vertex")
            elif w not in Q or vwLength <= Q[w]:
            	if w not in Q or vwLength < Q[w]:
                	Q[w] = vwLength
                	for index in range(len(list_P)):
                		list_P[index][w] = v
            	 # remove to disable ECMP
            	else:
            		if (ECMP == True):
            			ECMP_tuples1.append(w)
            			ECMP_tuples2.append(v)
    #print ("D = ", D)
    #print("list_P = ", list_P)
    #print(ECMP_tuples1, ECMP_tuples2)
    #print(len(ECMP_tuples))
    return (D,list_P, ECMP_tuples1, ECMP_tuples2)
    
                


def shortestPath(ECMP, G,start,end):
    list_Path = []
    D,list_P, ECMP_tuples1, ECMP_tuples2 = Dijkstra(ECMP, G,start,end)
    P = list_P[0]
    Path = []
    endP = end
    while 1:
       	Path.append(endP)
       	if endP == start: break
       	indices = [i for i in range(len(ECMP_tuples1)) if ECMP_tuples1[i] == endP]
       	#print(endP, indices)
       	if len(indices) > 0:
       		for index in indices:
       			P2 = P.copy()
       			P2[endP] = ECMP_tuples2[index]
       			list_PathN = path_loop(P2, Path.copy(), start, endP, ECMP_tuples1, ECMP_tuples2)
       			for new_path in list_PathN:
       				list_Path.append(new_path)
       	endP = P[endP]
    Path.reverse()
    list_Path.append(Path)
    
    #print(len(list_Path))
    return list_Path
    
def path_loop(P, Path, start, endP, ECMP_tuples1, ECMP_tuples2):
    list_Path = []
    while 1:
       	endP = P[endP]
       	Path.append(endP)
       	if endP == start: break
       	indices = [i for i in range(len(ECMP_tuples1)) if ECMP_tuples1[i] == endP]
       	#print(endP, indices)
       	if len(indices) > 0:
       		for index in indices:
       			P2 = P.copy()
       			P2[endP] = ECMP_tuples2[index]
       			list_PathN = path_loop(P2, Path.copy(), start, endP, ECMP_tuples1, ECMP_tuples2)
       			for new_path in list_PathN:
       				list_Path.append(new_path)
    Path.reverse()
    list_Path.append(Path)
    return list_Path
