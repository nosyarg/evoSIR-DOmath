import networkx as nx #this package is the backbone of the graph theory in this program
from random import *
from bisect import * #used to speed up array operations by taking advantage of sorting
def inlist(a, x):#determines an element x is in a sorted array a
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return 1
    return 0
while(1):#main loop
        #the next few lines determine the parameters for the graph on each run
        n = 100000
        mu = 5
        p = mu/n
        rho = 4#round(100*random())/10#round(100*random())/10
        lambd = 1.1 + .3*random()
        infected = []#this will store all the infected vertices
        si = []#this will store all the si edges
        G = nx.fast_gnp_random_graph(n,p)#set up the er graph and infect patient 0
        for i in range(len(G.nodes)):
                G.nodes[i]['state'] = 's'
        G.nodes[0]['state'] = 'i'
        infected = [0]
        si = sorted(list(G.edges(0)))#add all the edges on patient 0 to the si list
        '''
        This loop runs an exponential race at each timestep.
        It doesn't explicitly reference time, instead, it randomly selects which event will happen next.
        With probability i/(i+lambda*si+rho*si) that is an infected vertex becoming r.
        With probability (lambda*si)/(i+lambda*si+rho*si) that is an S becoming an I.
        With probability (rho*si)/(i+lambda*si+rho*si) this will be an si edge rewiring.
        '''
        while(len(infected) > 0): #This loop does evosir until there are no i's left
                if(random() < len(infected)/(len(infected)+(rho+lambd)*len(si))):#Recovery event
                        todie = choice(infected)#pick an infected vertex and kill it
                        G.nodes[todie]['state'] = 'r'
                        del infected[bisect_left(infected,todie)]
                        for node in G.neighbors(todie):#kill all the edges that used to be si connected to this vertex. Do it twice in case the edge is stored backwards
                                if(inlist(si,(todie,node))):
                                        del si[bisect_left(si,(todie,node))]
                                if(inlist(si,(node,todie))):
                                        del si[bisect_left(si,(node,todie))]
                elif(random()<lambd/(rho+lambd)): #infection edge event 
                        cross = choice(si)#pick the edge that the infection will cross
                        if(G.nodes[cross[0]]['state']=='s'):
                                G.nodes[cross[0]]['state'] = 'i'
                                insort(infected,cross[0])
                                for node in G.neighbors(cross[0]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[0],node))):
                                                        del si[bisect_left(si,(cross[0],node))]
                                                if(inlist(si,(node,cross[0]))):
                                                        del si[bisect_left(si,(node,cross[0]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[0],node))
                        elif(G.nodes[cross[1]]['state']=='s'):#do everything from the other side if the edge was reversed.
                                G.nodes[cross[1]]['state'] = 'i'
                                insort(infected,(cross[1]))
                                for node in G.neighbors(cross[1]):#change si list
                                        if(G.nodes[node]['state']=='i'):
                                                if(inlist(si,(cross[1],node))):
                                                        del si[bisect_left(si,(cross[1],node))]
                                                if(inlist(si,(node,cross[1]))):
                                                        del si[bisect_left(si,(node,cross[1]))]
                                        if(G.nodes[node]['state']=='s'):
                                                insort(si,(cross[1],node))
                else: #rewire event
                        rewire = choice(si)#grab the edge that you want to rewire
                        del si[bisect_left(si,rewire)]
                        G.remove_edge(rewire[0],rewire[1])#will this cause an issue? dont think so b/c 
                        if(G.nodes[rewire[0]]['state'] == 's'):#determine which of the vertices connected to the edge is susceptible
                                svert = rewire[0]
                        else:
                                svert = rewire[1]
                        newvert = svert
                        while((newvert==svert) or (newvert in G.neighbors(svert))):#this will loop until it selects a valid target for the rewire
                                newvert = int(n*random())
                        G.add_edge(newvert,svert)
                        if(G.nodes[newvert]['state'] == 'i'):
                                insort(si,(newvert,svert))
        survivors = 0
        for i in range(len(list(G.nodes))):#count the number of s vertices at the end
                if(G.nodes[i]['state']=='s'):
                        survivors += 1
        writefile = open('data/exptime.csv','a')#write the relevant data about the sim
        writefile.write(str(n) + ',' + str(p) + ',' + str(lambd) + ',' + str(rho) + ',' + str(survivors/n) + '\n')
        writefile.close()
        #print(survivors/n)
