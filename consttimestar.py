import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from random import *
from bisect import *
def inlist(a, x):
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return 1
    return 0
while(1):
        n = 100
        mu = 5
        lambd = 6
        rho = 4
        '''
        if(random() > .5):
                lambd = 1
                rho = 10*random()
        else:
                lambd = 10*random()
                rho = 4
                '''
        time = 0
        infected = []
        si = []
        G = nx.star_graph(n) #set up the er graph and infect patient 0
        for i in range(len(G.nodes)):
                G.nodes[i]['state'] = 's'
        G.nodes[0]['state'] = 'i'
        infected = [(0,time+1)]
        si = sorted(list(G.edges(0)))
        while(len(si) > 0): 
                timetoevent = np.random.exponential(1/(len(si)*(lambd+rho)))
                if(timetoevent+time>infected[0][1]):
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
                        cross = choice(si)
                        time += timetoevent
                        if(G.nodes[cross[0]]['state']=='s'):#the susceptible is first in the list
                                G.nodes[cross[0]]['state'] = 'i'
                                infected.append((cross[0],time+1))
                                for node in G.neighbors(cross[0]):#change si list
                                        if((G.nodes[node]['state']=='i')):#if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[0],node))):
                                                        del si[bisect_left(si,(cross[0],node))]
                                                if(inlist(si,(node,cross[0]))):
                                                        del si[bisect_left(si,(node,cross[0]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[0],node))
                        elif(G.nodes[cross[1]]['state']=='s'):#the susceptible is second in the list
                                G.nodes[cross[1]]['state'] = 'i'
                                infected.append((cross[1],time+1))
                                for node in G.neighbors(cross[1]):#change si list
                                        if((G.nodes[node]['state']=='i')):#if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[1],node))):
                                                        del si[bisect_left(si,(cross[1],node))]
                                                if(inlist(si,(node,cross[1]))):
                                                        del si[bisect_left(si,(node,cross[1]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[1],node))
                else: #rewire edge event
                        time += timetoevent
                        rewire = choice(si)
                        del si[bisect_left(si,rewire)]
                        if(G.nodes[rewire[0]]['state'] == 's'):
                                svert = rewire[0]
                        else:
                                svert = rewire[1]
                        newvert = svert
                        while((newvert==svert) or (newvert in G.neighbors(svert))):#pick a valid vertex
                                newvert = int(n*random())
                        G.remove_edge(rewire[0],rewire[1])
                        G.add_edge(newvert,svert)
                        if(G.nodes[newvert]['state'] == 'i'):
                                insort(si,(newvert,svert))
        survivors = 0
        for i in range(len(list(G.nodes))):
                if(G.nodes[i]['state']=='s'):
                        survivors += 1
        writefile = open('data/consttimestarcomponents.csv','a')
        writefile.write(str(n) + ',' + str(lambd) + ',' + str(rho) + ',' + str(survivors/n) + ',' + str(nx.number_connected_components(G)) + '\n')
        writefile.close()
        colorarr = [0]*(n+1)
        for i in range(n+1):
                if(G.nodes[i]['state'] == 's'):
                        colorarr[i] = 'green'
                elif(G.nodes[i]['state'] == 'i'):
                        colorarr[i] = 'blue'
                elif(G.nodes[i]['state'] == 'r'):
                        colorarr[i] = 'blue'
        nx.draw(G,node_size = 20, node_color = colorarr)
        plt.show()
        print(survivors/n)
        print(nx.number_connected_components(G))
