import networkx as nx
import numpy as np
from random import *
from bisect import *
def inlist(a, x):
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return 1
    return 0
n = 10000
survivors = n
while(survivors/n > .2):
        survivors = n
        writepath = 'data/consttimecritval.csv'
        writefile = open(writepath,'w')
        writefile.write('n,p,lambda,rho,s,i,r,mu,t\n')
        writefile.close()
        mu = 5
        susceptible = list(range(1,n))
        p = mu/n
        lambd = 1.0084#1#round(100*random())/10
        rho = 4#round(100*random())/10#4
        time = 0
        infected = []
        si = []
        G = 0
        G = nx.fast_gnp_random_graph(n,p) #set up the er graph and infect patient 0
        for i in range(len(G.nodes)):
                G.nodes[i]['state'] = 's'
        G.nodes[0]['state'] = 'i'
        infected = [(0,time+1)]
        si = sorted(list(G.edges(0)))
        counter = 0
        while(len(si) > 0): 
                counter += 1
                if(counter % 10 == 0):
                        mu = 0
                        for sus in susceptible:
                                mu += G.degree(sus)/len(susceptible)
                numinfected = len(infected)
                writefile = open(writepath,'a')
                writefile.write(str(n) + ',' + str(p) + ',' + str(lambd) + ',' + str(rho) + ',' + str(survivors/n) + ',' + str(numinfected/n) + ',' + str((n-numinfected - survivors)/n) +','+ str(mu) + ',' + str(time)  + '\n')
                writefile.close()
                timetoevent = np.random.exponential(1/(len(si)*(lambd+rho)))
                if(timetoevent+time>infected[0][1]):
                        #print('death')
                        survivors -= 1
                        time = infected[0][1]
                        todie = infected[0][0]
                        G.nodes[todie]['state'] = 'r'
                        del infected[0]#clean up infection lists
                        for node in G.neighbors(todie):
                                if(inlist(si,(todie,node))):
                                        del si[bisect_left(si,(todie,node))]
                                if(inlist(si,(node,todie))):
                                        del si[bisect_left(si,(node,todie))]
                elif(random()<lambd/(rho+lambd)): #infection edge event 
                        #print('infection')
                        cross = choice(si)
                        time += timetoevent
                        if(G.nodes[cross[0]]['state']=='s'):#the susceptible is first in the list
                                G.nodes[cross[0]]['state'] = 'i'
                                infected.append((cross[0],time+1))
                                del susceptible[bisect_left(susceptible,cross[0])]
                                for node in G.neighbors(cross[0]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[0],node))):
                                                        del si[bisect_left(si,(cross[0],node))]
                                                if(inlist(si,(node,cross[0]))):
                                                        del si[bisect_left(si,(node,cross[0]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[0],node))
                        elif(G.nodes[cross[1]]['state']=='s'):#the susceptible is second in the list
                                G.nodes[cross[1]]['state'] = 'i'
                                infected.append((cross[1],time+1))
                                del susceptible[bisect_left(susceptible,cross[1])]
                                for node in G.neighbors(cross[1]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[1],node))):
                                                        del si[bisect_left(si,(cross[1],node))]
                                                if(inlist(si,(node,cross[1]))):
                                                        del si[bisect_left(si,(node,cross[1]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[1],node))
                else: #rewire edge event
                        #print('rewire')
                        time += timetoevent
                        rewire = choice(si)
                        del si[bisect_left(si,rewire)]
                        G.remove_edge(rewire[0],rewire[1])
                        if(G.nodes[rewire[0]]['state'] == 's'):
                                svert = rewire[0]
                        else:
                                svert = rewire[1]
                        newvert = svert
                        while((newvert==svert) or (newvert in G.neighbors(svert))):#pick a valid vertex
                                newvert = int(n*random())
                        G.add_edge(newvert,svert)
                        if(G.nodes[newvert]['state'] == 'i'):
                                insort(si,(newvert,svert))
        print(survivors/n)
