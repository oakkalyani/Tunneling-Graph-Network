import zen
import sys
sys.path.append('../zend3js/')
#import d3js
import numpy
import time
from numpy import *
from time import sleep
import matplotlib.pyplot as plt
#import randomnet
import string
import math

print 'Start'

#viz = 0 #1 for visualizer


G = zen.io.gml.read('Plano2016.gml',weight_fxn = lambda x: x['weight'])
# zen.io.gml.write(G,'test.gml')
# print G.node_data[0]

n = G.num_nodes
A = G.matrix().transpose()
print n
print G.num_edges


'\n###########################################################\n'

print '\n\nCyclic / Acyclic Check'

def cycliccheck(G):
    A=G.matrix().transpose()
    OutM=A.sum(axis=0)
    n=len(A)
    cycheck=OutM.sum()
    if n!=0:
        for i in range(0,n):
            if OutM[i]==0:
                cycheck=-1
                break
        if cycheck==-1:
            G.rm_node_(i)
            G.compact()
            cycheck=cycliccheck(G)
    if cycheck==0:
        print '\nGraph is acyclic'
    else:
        print '\nGraph is cyclic'

cycle=cycliccheck(G)



print '\n###########################################################\n'
def print_top(G,v, num=5):
	idx_list = [(i,v[i]) for i in range(len(v))]
	idx_list = sorted(idx_list, key = lambda x: x[1], reverse=True)
	for i in range(min(num,len(idx_list))):
		nidx, score = idx_list[i]
		print '  %i. %s (%1.4f)' % (i+1,G.node_object(nidx),score)
		#print '  %i. %s' % (i+1,G.node_object(idx))

def index_of_max(v):
	return numpy.where(v == max(v))[0]

# Degree Centrality
print '\nDegree Centrality:'
vdeg=[0]*G.num_nodes
for i in range(1,n):
	vdeg[i]=sum(A[i])
print_top(G,vdeg)

print '\nEigenvector Centrality:'
k, V = numpy.linalg.eig(A)
k=numpy.abs(k)
k1_idx = index_of_max(k) # find the index of the largest eigenvalue
# finish printing the top 5 eigenvector centrality characters by linear algebra
vtemp=V[:,k1_idx]
veig=numpy.abs(vtemp)
veig=veig/sum(veig)
print_top(G,veig)

print '\nKatz Centrality:'
for i in range(1,6):
	a=1/float(i*i+max(k))	#to ensure alpha < 1/k1 condition
	b=numpy.ones((n,1))
	vkatz=numpy.dot(numpy.linalg.inv(numpy.eye(n)-a*A),b)
	print 'For alpha = %1.10f, Katz Centrality top 5 list is:' %a
	print_top(G,vkatz)

print '\nPageRank'
D=numpy.eye(n)
for i in range(0,n):
	D[i,i]=max(D[i,i],sum(A[:,i]))
kpr,vpr = numpy.linalg.eig(numpy.dot(A,numpy.linalg.inv(D)))
a=1/(1+max(kpr))
b=numpy.ones((n,1))
vpr=numpy.dot(numpy.dot(D,numpy.linalg.inv(D-a*A)),b)
# print_top(G,vpr)

print '\nBetweenness Centrality'
vbet=zen.algorithms.centrality.betweenness_centrality_(G,weighted=True)
print_top(G,vbet)

centrality = []
v = []
for i in range(0,n):
    # print (i)
    v.append(0.05*vdeg[i] + 0.4*veig[i] + 0.4*vbet[i] + 0.15*vkatz[i])
    v[i]=v[i][0]


print '\n###########################################################\n'

def degree_sequence(G):
	return [degree for degree,freq in enumerate(zen.degree.ddist(G,normalize=False)) for f in range(int(freq))]

def configuration_model(degree_sequence,G=None):
	import numpy.random as numpyrandom
	if G is None:
		G = zen.Graph()

	n = len(degree_sequence)
	for i in range(n):
		G.add_node(i)

	# this is twice the number of edges, needs to be even
	assert mod(sum(degree_sequence),2) == 0, 'The number of edges needs to be even; the sum of degrees is not even.'
	num_edges = sum(degree_sequence)/2

	# the number of edges should be even
	assert mod(num_edges,2) == 0, 'The number of edges needs to be even.'

	stubs = [nidx for nidx,degree in enumerate(degree_sequence) for d in range(degree)]
	stub_pairs = numpyrandom.permutation(num_edges*2)

	self_edges = 0
	multi_edges = 0
	for i in range(num_edges):
		uidx = stubs[stub_pairs[2*i]]
		vidx = stubs[stub_pairs[2*i+1]]
		if uidx == vidx:
			self_edges += 1
		if G.has_edge_(uidx,vidx):
			eidx = G.edge_idx_(uidx,vidx)
			G.set_weight_(eidx, G.weight_(eidx)+1 )
			multi_edges += 1
		else:
			G.add_edge_(uidx,vidx)

	print 'self edges: %i,  multi-edges: %i' % (self_edges,multi_edges)

	return G

def ddplot(ddist,k,ntype):
    plt.figure(figsize=(8,8))
    plt.suptitle('Degree Distribution: '+ntype)
    plt.bar(k,ddist, width=0.8, bottom=0, color='b')
    plt.xlabel('k')
    plt.ylabel('DegDist')
    plt.show()


#Given Network Properties
print 'Graph Network Fitting'

nn = G.num_nodes
ne = G.num_edges
print '\nBase model:'
print 'Number of Nodes: %i' %nn
print 'Number of Edges: %i' %ne
tp = 1

c = 2*ne/nn
Cg = zen.algorithms.clustering.gcc(G)

print 'Average degree: %i' %c
print 'Global Clustering coefficient: %1.4f' %Cg

ddist = zen.degree.ddist(G,normalize=False)
k = numpy.arange(len(ddist))

dseq = degree_sequence(G)
ntype='Base Model'
if tp==0:
    ddplot(ddist,k,ntype)

#Fitting ER
#print '\nErdos-Renyi Graph parameters to fit:'
#G2 = zen.Graph()
#pER = float(c)/(nn-1)
#print 'Edge existance probability: %1.5f' %pER
#G2 = randomnet.erdos_renyi(nn,pER,graph=G2)
#neER = G2.num_edges
#print 'Number of edges: %i' %neER
#cER = 2*G2.num_edges/nn
#CgER = zen.algorithms.clustering.gcc(G2)
#print 'Average degree: %i' %cER
#print 'Global Clustering coefficient: %1.4f' %CgER
#ddistER = zen.degree.ddist(G2,normalize=False)
#kER = numpy.arange(len(ddistER))
#ntype='ER'
#if tp==0:
#    ddplot(ddistER,kER,ntype)

#Fitting LA graph
print '\nLocal Attachment Graph parameters to fit:'
G2 = zen.DiGraph()
qLA = int(round(float(c)/2))
qrLA = 2
G2 = zen.generating.local_attachment(nn,qLA,qrLA,graph=G2)
neLA = G2.num_edges
print 'Number of outgoing edges from a node (q): %i' %qLA
print 'Number of initial outgoing edges from a node (qr): %i' %qrLA
print 'Number of edges: %i' %neLA
cLA = 2*G2.num_edges/nn
CgLA = zen.algorithms.clustering.gcc(G2)
print 'Average degree: %i' %cLA
print 'Global Clustering coefficient: %1.4f' %CgLA
ddistLA = zen.degree.ddist(G2,normalize=False)
kLA = numpy.arange(len(ddistLA))
ntype='LA'
if tp==0:
    ddplot(ddistLA,kLA,ntype)

#Fitting DD
print '\nDuplication Divergence Graph parameters to fit:'
G2 = zen.Graph()
pDD = 0.35
print 'Edge existance probability: %1.5f' %pDD
G2 = zen.generating.duplication_divergence_iky(nn,pDD,graph=G2)
neDD = G2.num_edges
print 'Number of edges: %i' %neDD
cDD = 2*G2.num_edges/nn
CgDD = zen.algorithms.clustering.gcc(G2)
print 'Average degree: %i' %cDD
print 'Global Clustering coefficient: %1.4f' %CgDD
ddistDD = zen.degree.ddist(G2,normalize=False)
kDD = numpy.arange(len(ddistDD))
ntype='DD'
if tp==0:
    ddplot(ddistDD,kDD,ntype)

#Fitting Config Model

print '\nConfiguration Model Graph parameters to fit:'
G2 = zen.Graph()
G2 = configuration_model(dseq,G=None)
neC = G2.num_edges
print 'Number of edges: %i' %neC
cC = 2*neC/nn
CgC = zen.algorithms.clustering.gcc(G2)
print 'Average degree: %i' %cC
print 'Global Clustering coefficient: %1.4f' %CgC
ddistC = zen.degree.ddist(G2,normalize=False)
kC = numpy.arange(len(ddistC))
ntype='Config Model'
if tp==0:
    ddplot(ddistC,kC,ntype)

print '\n###########################################################\n'
#Modularity and Communities
print 'Identifying node Communities in the Graph'
def modularity(G,c):
	d = dict()
	for k,v in c.iteritems():
		for n in v:
			d[n] = k
	Q, Qmax = 0,1
	for u in G.nodes_iter():
		for v in G.nodes_iter():
			if d[u] == d[v]:
				Q += ( int(G.has_edge(v,u)) - G.in_degree(u)*G.out_degree(v)/float(G.num_edges) )/float(G.num_edges)
				Qmax -= ( G.in_degree(u)*G.out_degree(v)/float(G.num_edges) )/float(G.num_edges)
	return Q, Qmax

def scalar_assortativity(G,d):
    x = zeros(G.num_nodes)
    for i in range(G.num_nodes):
        x[i] = d[G.node_object(i)]
    return scalar_assortativity_(G,x)

def scalar_assortativity_(G,x):
    A = G.matrix().T
    M = 2*A.sum().sum()
    ki = A.sum(axis=1) #row sum is in-degree
    ko = A.sum(axis=0) #column sum is out-degree
    mu = ( dot(ki,x)+dot(ko,x) )/M
    print 'Mean: %1.4f' %mu

    R, Rmax = 0, 0
    for i in range(G.num_nodes):
        for j in range(G.num_nodes):
             R += ( A[i,j]*(x[i]-mu)*(x[j]-mu) )/M
             Rmax += ( A[i,j]*(x[i]-mu)**2 )/M

    return R, Rmax

cset = zen.algorithms.community.louvain(G)
cnum = cset.__len__()
csize = zeros(cnum)

f = open('CommSet.txt','w')
for i in range(0,cnum):
    c = cset.community(i)
    csize[i] = c.__len__()
    cn = c.nodes()
    for j in range(0,len(cn)):
        f.write(str(cn[j]+1)+' ')
    f.write('\n')
f.close()

print '\n###########################################################\n'
## Visualizer
#def visualizer1(G):
#    x0=0
#    y0=0
#    d3 = d3js.D3jsRenderer(G, canvas_size=(2000,2000),interactive=False, autolaunch=False)
#    scale = float(2000)/60000
#    for i in range(0,G.num_nodes):
#        coord = G.node_data_(i)['zenData']
#        coord['xcord'] = coord['xcord']*scale
#        coord['ycord'] = coord['ycord']*scale
#        print coord
#        d3.position_node_(i,coord['xcord']+x0,coord['ycord']+y0)
#    d3.update()
#    d3.stop_server()
#
#if viz==1:
#    visualizer1(G)
#
#'\n###########################################################\n'
#def visualizer2(G,path):
#    d3 = d3js.D3jsRenderer(G, canvas_size=(2000,2000),interactive=False, autolaunch=False)
#    x0=0
#    y0=0
#    scale = float(2000)/60000
#    for i in range(0,G.num_nodes):
#        coord = G.node_data_(i)['zenData']
#        coord['xcord'] = coord['xcord']*scale
#        coord['ycord'] = coord['ycord']*scale
#        # print coord
#        d3.position_node_(i,coord['xcord']+x0,coord['ycord']+y0)
#    d3.highlight_edges_(path)
#    d3.update()
#    d3.stop_server()
#
#def visualizer3(G,path):
#    d3 = d3js.D3jsRenderer(G, canvas_size=(2000,2000),interactive=False, autolaunch=False)
#    x0=0
#    y0=0
#    scale = float(2000)/60000
#    for i in range(0,G.num_nodes):
#        coord = G.node_data_(i)['zenData']['zenData']
#        coord['xcord'] = coord['xcord']*scale
#        coord['ycord'] = coord['ycord']*scale
#        # print coord
#        d3.position_node_(i,coord['xcord']+x0,coord['ycord']+y0)
#    d3.highlight_edges_(path)
#    d3.update()
#    d3.stop_server()
##Diameter
# d = zen.diameter(G)
print 'Shortest Path and Cut-Set analysis for given road map'
D,P = zen.algorithms.shortest_path.all_pairs_dijkstra_path_(G)
m=D[0][0]
uobj_w = []
vobj_w = []
for i in range(0,G.num_nodes):
    for j in range(i+1,G.num_nodes):
        if m<D[i][j]:
            m=D[i][j]
            uobj_w = []
            vobj_w = []
            uobj_w.append(G.node_object(i))
            vobj_w.append(G.node_object(j))
        elif m==D[i][j]:
            uobj_w.append(G.node_object(i))
            vobj_w.append(G.node_object(j))
print 'Maximum of weighted shortest paths: %i : Between %i and %i' %(m,uobj_w[0],vobj_w[0])

path_w = zen.algorithms.shortest_path.pred2path_(G.node_idx(uobj_w[0]),G.node_idx(vobj_w[0]),P)
print 'Node Path:'
for i in range(0,len(path_w)):
    path_w[i]=G.node_object(path_w[i])
    print path_w[i],
    if i!=len(path_w)-1:
        print ' - ',
print ' '

vpath = []
for i in range(0,len(path_w)-1):
    vpath.append(G.edge_idx(path_w[i],path_w[i+1]))
if viz==1:
    visualizer2(G,vpath)
    time.sleep(3)

D,P = zen.algorithms.shortest_path.all_pairs_dijkstra_path_(G,ignore_weights=True)
m=D[0][0]
uobj_uw = []
vobj_uw = []
for i in range(0,G.num_nodes):
    for j in range(i+1,G.num_nodes):
        if m<D[i][j] and D[i][j]>0 and math.isinf(D[i][j])==0:
            m=D[i][j]
            uobj_uw = []
            vobj_uw = []
            uobj_uw.append(G.node_object(i))
            vobj_uw.append(G.node_object(j))
        elif m==D[i][j] and D[i][j]>0 and math.isinf(D[i][j])==0:
            uobj_uw.append(G.node_object(i))
            vobj_uw.append(G.node_object(j))
print 'Maximum of unweighted shortest paths: %i : Between %i and %i' %(m,uobj_uw[0],vobj_uw[0])

path_uw = zen.algorithms.shortest_path.pred2path_(G.node_idx(uobj_uw[0]),G.node_idx(vobj_uw[0]),P)
print 'Node Path:'
for i in range(0,len(path_uw)):
    path_uw[i]=G.node_object(path_uw[i])
    print path_uw[i],
    if i!=len(path_uw)-1:
        print ' - ',
print ' '

vpath = []
for i in range(0,len(path_uw)-1):
    vpath.append(G.edge_idx(path_uw[i],path_uw[i+1]))
if viz==1:
    visualizer2(G,vpath)
    time.sleep(3)


G2 = zen.DiGraph()
for i in range (0,G.num_nodes):
    G2.add_node(G.node_object(i),G.node_data_(i))
for i in range (0,G.num_nodes):
    # if G.is_valid_node_idx(i)!=1:
    #     print '%i is not valid' %i
    # else:
    #     print '%i is valid' %i
    nbrs = G.neighbors_(i)
    if len(nbrs)>0:
        for j in range(i+1,G.num_nodes):
            n1 = G.node_object(i)
            n2 = G.node_object(j)
            # print '%i - %i' %(n1,n2)
            if G.has_edge(n1,n2)==1:
                w = G.weight(n1,n2)
                # print w
                G2.add_edge(n1,n2,weight=w)
                G2.add_edge(n2,n1,weight=w)
                # print G2.weight_(G2.edge_idx(n1,n2))

path_cutw =zen.algorithms.flow.min_cut_set_(G2,G2.node_idx(uobj_w[0]),G2.node_idx(vobj_w[0]), capacity='weight')
print 'Weighted network cutset path(s):'
for i in range(0,len(path_cutw)):
    a,b = G2.endpoints(path_cutw[i])
    print '     %i - %i' %(a,b)
# a=G2.node_object(a)
# b=G2.node_object(b)
# print 'Path: %i - %i' %(G2.node_idx(uobj_w[0]),G2.node_idx(vobj_w[0]))
# print 'Cut-set: %i - %i' %(a,b)

path_cutuw =zen.algorithms.flow.min_cut_set_(G2,G2.node_idx(uobj_uw[0]),G2.node_idx(vobj_uw[0]), capacity='unit')
print 'Unweighted network cutset path(s):'
for i in range(0,len(path_cutuw)):
    a,b = G2.endpoints(path_cutuw[i])
    print '     %i - %i' %(a,b)
# a,b = G2.endpoints_(path_cutuw[0])
# a=G2.node_object(a)
# b=G2.node_object(b)
# print 'Path: %i - %i' %(G2.node_idx(uobj_uw[0]),G2.node_idx(vobj_uw[0]))
# print 'Cut-set: %i - %i' %(a,b)

print '\n###########################################################\n'
########################################################################################################
##Analysis for main body of nodes removing outer hanging edges

print 'Shortest path and Cut-set analysis for graph with no hanging edges'

G3 = G
# print G3.node_data_(0)
# for i in range(0,G3.num_nodes):
#     G3.set_node_data(i,G3.node_data_(i)['zenData'])
# print G3.node_data_(i)
check = 1
while check!=0:
    check = 0
    for i in range(0,G3.num_nodes):
        if G3.degree_(i)==1:
            G3.rm_edge_(G3.edge_idx_(G3.neighbors_(i)[0],i))
            check=check+1
time.sleep(3)
if viz==1:
    visualizer1(G3)

D,P = zen.algorithms.shortest_path.all_pairs_dijkstra_path_(G3)
# print D
m=0
uobj_w = []
vobj_w = []
for i in range(0,G3.num_nodes):
    for j in range(i+1,G3.num_nodes):
        if m<D[i][j] and D[i][j]>0 and math.isinf(D[i][j])==0:
            m=D[i][j]
            uobj_w = []
            vobj_w = []
            uobj_w.append(G3.node_object(i))
            vobj_w.append(G3.node_object(j))
        elif m==D[i][j] and D[i][j]>0 and math.isinf(D[i][j])==0:
            uobj_w.append(G3.node_object(i))
            vobj_w.append(G3.node_object(j))
# print uobj_w[0]
# print vobj_w[0]
print 'Maximum of weighted shortest paths: %i : Between %i and %i' %(m,uobj_w[0],vobj_w[0])
# print m
# print uobj_w
# print vobj_w
path_w = zen.algorithms.shortest_path.pred2path_(G3.node_idx(uobj_w[0]),G3.node_idx(vobj_w[0]),P)
p = len(path_w)
for i in range(0,p):
    path_w[i]=G3.node_object(path_w[i])
    print path_w[i],
    if i!=p-1:
        print ' - ',
print ' '

vpath = []
for i in range(0,p-1):
    vpath.append(G3.edge_idx(path_w[i],path_w[i+1]))

if viz==1:
    visualizer2(G3,vpath)
    time.sleep(3)
# print path_uw
# print D
# print P
# print d
# print uidx
# print vidx
D,P = zen.algorithms.shortest_path.all_pairs_dijkstra_path_(G3,ignore_weights=True)
m=D[0][0]
uobj_uw = []
vobj_uw = []
for i in range(0,G3.num_nodes):
    for j in range(i+1,G3.num_nodes):
        if m<D[i][j] and D[i][j]>0 and math.isinf(D[i][j])==0:
            m=D[i][j]
            uobj_uw = []
            vobj_uw = []
            uobj_uw.append(G3.node_object(i))
            vobj_uw.append(G3.node_object(j))
        elif m==D[i][j] and D[i][j]>0 and math.isinf(D[i][j])==0:
            uobj_uw.append(G3.node_object(i))
            vobj_uw.append(G3.node_object(j))
print 'Maximum of unweighted shortest paths: %i : Between %i and %i' %(m,uobj_uw[0],vobj_uw[0])
# print m
# print uobj_uw
# print vobj_uw
path_uw = zen.algorithms.shortest_path.pred2path_(G3.node_idx(uobj_uw[0]),G3.node_idx(vobj_uw[0]),P)
p = len(path_uw)
for i in range(0,p):
    path_uw[i]=G3.node_object(path_uw[i])
    print path_uw[i],
    if i!=p-1:
        print ' - ',
print ' '
# print path_uw
vpath = []
for i in range(0,p-1):
    vpath.append(G3.edge_idx(path_uw[i],path_uw[i+1]))

if viz==1:
    visualizer2(G3,vpath)
    time.sleep(3)

G2 = zen.DiGraph()
for i in range (0,G3.num_nodes):
    G2.add_node(G3.node_object(i),G3.node_data_(i))
for i in range (0,G3.num_nodes):
    # if G3.is_valid_node_idx(i)!=1:
    #     print '%i is not valid' %i
    # else:
    #     print '%i is valid' %i
    nbrs = G3.neighbors_(i)
    if len(nbrs)>0:
        for j in range(i+1,G3.num_nodes):
            n1 = G3.node_object(i)
            n2 = G3.node_object(j)
            if G3.has_edge(n1,n2)==1:
                w = G3.weight(n1,n2)
                # print w
                G2.add_edge(n1,n2,weight=w)
                G2.add_edge(n2,n1,weight=w)
                # print G2.weight_(G2.edge_idx(n1,n2))

# print uobj_w[0]
# print G2.node_idx(uobj_w[0])
# print G2.node_idx(vobj_w[0])
path_cutw =zen.algorithms.flow.min_cut_set_(G2,G2.node_idx(uobj_w[0]),G2.node_idx(vobj_w[0]), capacity='weight')
print 'Weighted network cutset path(s):'
for i in range(0,len(path_cutw)):
    a,b = G2.endpoints(path_cutw[i])
    print '     %i - %i' %(a,b)

path_cutuw =zen.algorithms.flow.min_cut_set_(G2,G2.node_idx(uobj_uw[0]),G2.node_idx(vobj_uw[0]), capacity='unit')
print 'Unweighted network cutset path(s):'
for i in range(0,len(path_cutuw)):
    a,b = G2.endpoints(path_cutuw[i])
    print '     %i - %i' %(a,b)


print '\n###########################################################\n'

#Bridge Map generation
cen_list = []
for i in range(0,G.num_nodes):
    cen_list.append([G.node_object(i),v[i]])
print cen_list[0][0]
print cen_list[0][1]
cen_rank=sorted(cen_list,key=lambda x:x[1],reverse=True)

# print cen_rank
print len(cen_rank)
print len(cen_rank[0])
vin = 10
Gb = zen.Graph()
for i in range(0,vin):
    Gb.add_node()

print Gb.num_nodes
