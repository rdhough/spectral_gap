import numpy as np
from scipy import signal
from scipy import misc
import bisect
from copy import copy, deepcopy
from scipy.optimize import linprog
import itertools as it

target_value = 5.9776
target_value_1 = 5.9776- 1.35
target_value_2 = 5.9776 - 2.7
target_value_3 = 5.9776 - 3.5
target_value_4 = 5.9776 - 4

# Given a list of configurations with all nodes of height 1, builds all configurations that can be obtained by adding one node.

def GetCandidates(l):
    Candidates = []
    for c in l:
        for n in c.NextGeneration():
            if NotPresent(Candidates, n):
                bisect.insort(Candidates, n)
    return Candidates

# Given a list of configurations with a node of height 2, builds all configurations that can be obtained by adding one node

def GetCandidates2(l):
    Candidates = []
    for c in l:
        for n in c.NextGeneration2():
            if NotPresent(Candidates, n):
                bisect.insort(Candidates, n)
    return Candidates

# Given a list of configurations with two adjacent nodes of height 2, builds all configurations that can be obtained by adding one node

def GetCandidates22(l):
    Candidates = []
    for c in l:
        for n in c.NextGeneration22():
            if NotPresent(Candidates, n):
                bisect.insort(Candidates, n)
    return Candidates

# In the recurssion building components of height 1, given a list of configurations of a certain size, generate the next size configurations

def GenNextLevel(l, numbins):
    counter = 0;
    CurrentLevel = []
    NextLevel = []
    for c in l:
        counter +=1
        min_value = 10
        for cc in c.NextGeneration():
            cc.NumBins = numbins
            v = cc.ObtainValueCareful1()
            if v < target_value:
                if NotPresent(NextLevel, cc):
                    bisect.insort(NextLevel, cc)
            min_value = min(min_value, v)
        if min_value < target_value:
            bisect.insort(CurrentLevel, c)
    return [CurrentLevel, NextLevel]

# Recursively build all admissible components with nodes of height 1

def RecursivelyBuildComponents():
    v = Vertex(0,0,0)
    c = Configuration([v])
    c.ObtainOptimizationValue()
    ConfigurationList = []
    ConfigurationList.append([c])
    counter = 1
    while len(ConfigurationList[-1]) > 0 and counter < 6:
        counter+=1
        candidates = GetCandidates(ConfigurationList[-1])
        newlist = []
        for c in candidates:
            ov = c.ObtainValueCareful1()
            if ov < target_value:
                bisect.insort(newlist, c)
        ConfigurationList.append(newlist)
    return ConfigurationList

# Recursively build all admissable components with a node of height 2

def RecursivelyBuildComponents2():
    v0 = Vertex(0,0,0)
    c = Configuration([v0])
    c.LaplaceConstraint = [2]
    l = [c]
    ConfigurationList = [l]
    while len(ConfigurationList[-1]) > 0:
        candidates = GetCandidates2(ConfigurationList[-1])
        newlist = []
        for c in candidates:
            if NotPresent(newlist, c):
                ov = c.ObtainValueCareful2()
                if ov < target_value:
                    bisect.insort(newlist, c)
        ConfigurationList.append(newlist)
    return ConfigurationList

# Recursively build all admissable components with two adjacent nodes of height 2 

def RecursivelyBuildComponents22():
    v0 = Vertex(0,0,0)
    v1 = Vertex(0,0,1)
    c = Configuration([v0,v1])
    c.LaplaceConstraint = [2,2]
    l = [c]
    ConfigurationList = [l]
    while len(ConfigurationList[-1]) > 0:
        candidates = GetCandidates22(ConfigurationList[-1])
        newlist = []
        for c in candidates:
            if NotPresent(newlist, c):
                ov = c.ObtainValueCareful2()
                if ov < target_value:
                    bisect.insort(newlist, c)
        ConfigurationList.append(newlist)
    return ConfigurationList
    
# Given a list of configurations of size 4, tests all assignments of signs such that the resulting configuration is C2

def Obtain4C2Components(l):
    returnlist = []
    for c in l:
        vx = 0
        vy = 0
        offset = 0
        for v in c.vertices:
            vx += v[1]*.5
            vy += v[2]*.5
            offset += v[3]*.5
        v_ = Vertex(vx - c.vertices[0][1], vy-c.vertices[0][2], offset-c.vertices[0][3])
        if not NotPresent(c.vertices, v_):
            ind = bisect.bisect(c.vertices, v_) - 1
            laplace = [-1,1,1,1]
            laplace[ind] = -1
            cc = Configuration(c.vertices)
            cc.setLaplacian(laplace)
            cc.setLaplaceConstraint(laplace)
            returnlist.append(cc)
    return returnlist

# Given a list of configurations of size 6, tests all assignments of signs such that the resulting configuration is C2

def Obtain6C2Components(l):
    returnlist = []
    for c in l:
        vx = 0
        vy = 0
        offset = 0
        for v in c.vertices:
            vx += 3*v[1] + v[3]                 
            vy += 3*v[2] + v[3]                 
        verts = c.vertices
        for i in range(1,5):
            for j in range(i+1,6):
                if (6*(verts[0][1] + verts[i][1] + verts[j][1])+2*(verts[0][3] + verts[i][3] + verts[j][3]) == vx) and (6*(verts[0][2] + verts[i][2] + verts[j][2])+2*(verts[0][3] + verts[i][3] + verts[j][3]) == vy) :
                    laplace=[-1,1,1,1,1,1]
                    laplace[i] = -1
                    laplace[j] = -1
                    cc = Configuration(c.vertices)
                    cc.setLaplacian(laplace)
                    cc.setLaplaceConstraint(laplace)
                    returnlist.append(cc)
    return returnlist

# Given a list of configurations of size 8, tests all assignments of signs such that the resulting configuration is C2

def Obtain8C2Components(l):
    returnlist = []
    for c in l:
        vx = 0
        vy = 0
        offset = 0
        for v in c.vertices:
            vx += 3*v[1]+ v[3] 
            vy += 3*v[2] + v[3] 
        print c
        print [vx,vy]
        verts = c.vertices
        for i in range(1,6):
            for j in range(i+1,7):
                for k in range(j+1,8):
                    if (6*(verts[0][1] + verts[i][1] + verts[j][1]+verts[k][1])+2*(verts[0][3] + verts[i][3] + verts[j][3] + verts[k][3]) == vx) and (6*(verts[0][2] + verts[i][2] + verts[j][2]+verts[k][2])+2*(verts[0][3] + verts[i][3] + verts[j][3]+verts[k][3]) == vy) :
                        laplace=[-1,1,1,1,1,1,1,1]
                        laplace[i] = -1
                        laplace[j] = -1
                        laplace[k] = -1
                        cc = Configuration(c.vertices)
                        cc.setLaplacian(laplace)
                        cc.setLaplaceConstraint(laplace)
                        returnlist.append(cc)
    return returnlist

# Given a list of configurations, iterates over the list, and all nodes which could be assigned height 2, and tests the resulting value

def ObtainHeight2Configs(l):
    return_list = []
    for c in l:
        for j in range(c.n):
            c.LaplaceConstraint = np.ones(c.n)
            c.LaplaceConstraint[j] = 2
            val = c.ObtainValueCareful2()
            if val < target_value:
                return_list.append([c, c.LaplaceConstraint, val])
    return return_list

# Given a configuration of size 3, iterates over all assignments of signs, and determines the location of a singleton such that the configuration is in C2

def FindSize3SingletonPair(c):
    returnlist = []
    for i_ in range(2):
        c1 = Configuration(c.vertices)
        c1.Laplacian = np.ones(3)
        c1.Laplacian[i_] = -1
        v1 = c1.vertices[i_]
        vx = 0
        vy = 0
        vz = 0
        for m in range(3):
            vx += c1.Laplacian[m] * c1.vertices[m].x_coor
            vy += c1.Laplacian[m] * c1.vertices[m].y_coor
            vz += c1.Laplacian[m] * c1.vertices[m].parity
        if vz % 3 != 2:
            if vz %3 == 0:
                vx += vz//3
                vy += vz//3
                vz = 0
            elif vz%3 == 1:
                vx += (vz-1)//3
                vy += (vz-1)//3
                vz = 1
            v = Vertex(int(vx), int(vy), int(vz))
            if NotPresent(c1.vertices, v) and NotPresent(c1.dist2neighbors, v):
                c2 = c1.AddElem(v)
                c2.Laplacian = np.ones(4)
                j0 = bisect.bisect_left(c2.vertices, v)
                j1 = bisect.bisect_left(c2.vertices, v1)
                c2.Laplacian[j0] = -1
                c2.Laplacian[j1] = -1
                returnlist.append(c2)
    return returnlist
                
# Given a configuration of size 5, iterates over all assignments of signs and determines the corresponding location of a singleton such that the configuration is in C2

def FindSize5SingletonPair(c):
    returnlist = []
    for minus_indices in it.combinations(range(4),2):
        c1 = c
        c1.Laplacian = np.ones(5)
        c1.Laplacian[minus_indices[0]] = -1
        c1.Laplacian[minus_indices[1]] = -1
        v1 = c1.vertices[minus_indices[0]]
        v2 = c1.vertices[minus_indices[1]]
        vx = 0
        vy = 0
        vz = 0
        for m in range(5):
            vx += c1.Laplacian[m] * c1.vertices[m].x_coor
            vy += c1.Laplacian[m] * c1.vertices[m].y_coor
            vz += c1.Laplacian[m] * c1.vertices[m].parity
        if vz % 3 != 2:
            if vz %3 == 0:
                vx += vz//3
                vy += vz//3
                vz = 0
            elif vz %3 == 1:
                vx += (vz-1)//3
                vy += (vz-1)//3
                vz = 1
            v = Vertex(int(vx), int(vy), int(vz))
            if NotPresent(c1.vertices, v) and NotPresent(c1.dist2neighbors, v):
                c2 = c1.AddElem(v)
                c2.Laplacian = np.ones(6)
                i0 = bisect.bisect_left(c2.vertices, v)
                i1 = bisect.bisect_left(c2.vertices, v1)
                i2 = bisect.bisect_left(c2.vertices, v2)
                c2.Laplacian[i0]=-1
                c2.Laplacian[i1]=-1
                c2.Laplacian[i2]=-1
                c2.LaplaceConstraint[i0]=-1
                c2.LaplaceConstraint[i1]=-1
                c2.LaplaceConstraint[i2]=-1
                bisect.insort(returnlist, c2)
    return returnlist

# This function tests whether the Laplacian of a configuration makes it mean 0

def MeanZero(c):
    laplacian = c.Laplacian
    vertices = c.vertices
    vx = 0
    vy = 0
    vz = 0
    for i in range(len(vertices)):
        vx += laplacian[i] * vertices[i].x_coor
        vy += laplacian[i] * vertices[i].y_coor
        vz += laplacian[i] * vertices[i].parity
    return (3*vx + vz == 0) and (3*vy + vz == 0)

# Given a list of configurations of size 5, finds all assignments of Laplacian values such that 0 has height 2, any adjacent nodes to 0 have value -1, and such that the configuration is in C2

def Obtain5C2_2(l):
    returnlist = []
    v0 = Vertex(0,0,0)
    v1 = Vertex(0,0,1)
    v2 = Vertex(-1,0,1)
    v3 = Vertex(0,-1,1)
    v_list = [v0, v1, v2, v3]
    for c in l:
        indices = []
        for v in v_list:
            if not NotPresent(c.vertices, v):
                i = bisect.bisect_left(c.vertices, v)
                bisect.insort(indices, i)
        i = bisect.bisect_left(c.vertices, v0)
        extra_indices = []
        for k in range(5):
            if NotPresent(indices, k):
                extra_indices.append(k)
        num_neg = 4 - len(indices)
        for neg_indices in it.combinations(extra_indices, num_neg):
            c0 = c
            c0.LaplaceConstraint = np.ones(5)
            c0.Laplacian = np.ones(5)
            c0.LaplaceConstraint[i] = 2
            c0.Laplacian[i] = 2
            for j in indices:
                if j != i:
                    c0.LaplaceConstraint[j] = -1
                    c0.Laplacian[j] = -1
            for j in neg_indices:
                c0.LaplaceConstraint[j] = -1
                c0.Laplacian[j] = -1
            print c0.Laplacian
            if MeanZero(c0):
                returnlist.append(c0)
    return returnlist

# Given a list of configurations of size 7, finds all assignments of Laplacian values such that 0 has height 2, any adjacent nodes to 0 have value -1, and such that the configuration is in C2

def Obtain7C2_2(l):
    returnlist = []
    v0 = Vertex(0,0,0)
    v1 = Vertex(0,0,1)
    v2 = Vertex(-1,0,1)
    v3 = Vertex(0,-1,1)
    v_list = [v0, v1, v2, v3]
    for c in l:
        indices = []
        for v in v_list:
            if not NotPresent(c.vertices, v):
                i = bisect.bisect_left(c.vertices, v)
                bisect.insort(indices, i)
        i = bisect.bisect_left(c.vertices, v0)
        extra_indices = []
        for k in range(7):
            if NotPresent(indices, k):
                extra_indices.append(k)
        num_neg = 5 - len(indices)
        for neg_indices in it.combinations(extra_indices, num_neg):
            c0 = Configuration(c.vertices)
            c0.LaplaceConstraint = np.ones(7)
            c0.Laplacian = np.ones(7)
            c0.LaplaceConstraint[i] = 2
            c0.Laplacian[i] = 2
            for j in indices:
                if j != i:
                    c0.LaplaceConstraint[j] = -1
                    c0.Laplacian[j] = -1
            for j in neg_indices:
                c0.LaplaceConstraint[j] = -1
                c0.Laplacian[j] = -1
            print c0.Laplacian
            if MeanZero(c0):
                returnlist.append(c0)
    return returnlist

# Given a list of configurations of size 4 with nodes of height 2 at 0 and v, finds all ways of assigning signs so that the configuration is in C2

def Obtain4C2_22(l):
    return_list = []
    v0 = Vertex(0,0,0)
    v1 = Vertex(0,0,1)
    for c in l:
        i = bisect.bisect_left(c.vertices, v0)
        j = bisect.bisect_left(c.vertices, v1)
        indices = []
        for k in range(4):
            if k != i and k != j:
                indices.append(k)
        for count in range(2):
            if count == 0:
                c1 = c
                c1.LaplaceConstraint[indices[0]] = 1
                c1.LaplaceConstraint[indices[1]] = -1
                c1.LaplaceConstraint[i] = 2
                c1.LaplaceConstraint[j] = -2
                c1.Laplacian[indices[0]] = 1
                c1.Laplacian[indices[1]] = -1
                c1.Laplacian[i]=2
                c1.Laplacian[j]=-2
                if MeanZero(c1):
                    return_list.append(c1)
            if count == 1:
                c2 = c
                c2.LaplaceConstraint[indices[0]]=-1
                c2.LaplaceConstraint[indices[1]]=1
                c2.LaplaceConstraint[i]=2
                c2.LaplaceConstraint[j]=-2
                c2.Laplacian[indices[0]]=-1
                c2.Laplacian[indices[1]]=1
                c2.Laplacian[i]=2
                c2.Laplacian[j]=-2
                if MeanZero(c2):
                    return_list.append(c2)
    return return_list

# Given a list of configurations of size 6 with nodes at 0, v of height 2, this function finds all ways of assigning signs to the vertices so that the configuration is in C2

def Obtain6C2_22(l):
    return_list = []
    v0 = Vertex(0,0,0)
    v1 = Vertex(0,0,1)
    for c in l:
        i = bisect.bisect_left(c.vertices, v0)
        j = bisect.bisect_left(c.vertices, v1)
        indices = []
        for k in range(6):
            if k != i and k != j:
                indices.append(k)
        for neg_vars in it.combinations(range(4),2):
            c1 = c
            c1.LaplaceConstraint = np.ones(6)
            c1.Laplacian = np.ones(6)
            c1.LaplaceConstraint[i] = 2
            c1.LaplaceConstraint[j] = -2
            for m in neg_vars:
                c1.LaplaceConstraint[indices[m]] = -1
            c1.Laplacian[i] = 2
            c1.Laplacian[j] = -2
            for m in neg_vars:
                c1.Laplacian[indices[m]] = -1
            if MeanZero(c1):
                return_list.append(c1)
    return return_list

# This function enumerates all pairs of size two configurations both of which are nodes at distance 2 from each other.  The pairs are required to have the same orientation.

def findDist2Pairs():
    return_list = []
    v0 = Vertex(0,0,0)
    v1 = Vertex(1,0,0)
    for v2 in v0.ObtainDist3Neighbors():
        v3 = Vertex(v2.x_coor+1, v2.y_coor, v2.parity)
        v4 = Vertex(v2.x_coor-1, v2.y_coor, v2.parity)
        if NotPresent(v0.ObtainDist2Neighbors(), v2) and NotPresent(v1.ObtainDist2Neighbors(), v2):
            if NotPresent(v0.ObtainDist2Neighbors(), v3) and NotPresent(v1.ObtainDist2Neighbors(), v3):
                c1 = Configuration([v0,v1,v2,v3])
                if NotPresent(return_list, c1):
                    bisect.insort(return_list, c1)
            if NotPresent(v0.ObtainDist2Neighbors(), v4) and NotPresent(v1.ObtainDist2Neighbors(), v4):
                c2 = Configuration([v0,v1,v2,v4])
                if NotPresent(return_list, c2):
                    bisect.insort(return_list, c2)
    return return_list

# This function enumerates all sets of two adjacent vertices at distance 3 from each other

def TwoPairsAdjacent():
    return_list = []
    v0 = Vertex(0,0,0)
    v1 = Vertex(0,0,1)
    for v2 in v0.ObtainDist3Neighbors():
        if NotPresent(v0.ObtainDist2Neighbors(), v2) and NotPresent(v1.ObtainDist2Neighbors(), v2):
            for v3 in v2.ObtainNeighbors():
                if NotPresent(v0.ObtainDist2Neighbors(), v3) and NotPresent(v1.ObtainDist2Neighbors(), v3):
                    c = Configuration([v0,v1,v2,v3])
                    if NotPresent(return_list, c):
                        bisect.insort(return_list, c)
    return return_list


# This function enumerates all pairs of size two configurations, each of which is a pair of adjacent nodes.  The pairs are required to have the same orientation.

def findAdjacentPairs():
    return_list = []
    v0 = Vertex(0,0,0)
    v1 = Vertex(0,0,1)
    for v2 in v0.ObtainDist2TriangularLatticeNeighbors():
        v3 = Vertex(v2.x_coor, v2.y_coor, 1)
        c = Configuration([v0,v1,v2,v3])
        bisect.insort(return_list, c)
    return return_list

# Given a list of size 4 configurations this function finds all ways of attaching a pair of vertices at distance 3

def size4adjacentpairs(l):
    returnlist = []
    for c in l:
        candidate_list = []
        for v in c.vertices:
            candidate_list = MergeList(candidate_list, v.ObtainDist3Neighbors())
        for v1 in candidate_list:
            if NotPresent(c.vertices, v1) and NotPresent(c.dist2neighbors, v1):
                for v2 in v1.ObtainNeighbors():
                    if NotPresent(c.vertices, v2) and NotPresent(c.dist2neighbors, v2):
                        c1 = c.AddElem(v1)
                        c2 = c1.AddElem(v2)
                        bisect.insort(returnlist, c2)
    return returnlist

# Given a list of size 4 configurations this function finds all ways of attaching a pair of vertices which are 2-spaced

def size4dist2pairs(l):
    returnlist = []
    for c in l:
        candidate_list = []
        for v in c.vertices:
            candidate_list = MergeList(candidate_list, v.ObtainDist3Neighbors())
        for v1 in candidate_list:
            if NotPresent(c.vertices, v1) and NotPresent(c.dist2neighbors,v1):
                for v2 in v1.ObtainDist2Neighbors():
                    if NotPresent(c.vertices, v2) and NotPresent(c.dist2neighbors, v2) and NotPresent(v1.ObtainNeighbors(), v2) and v1 != v2:
                        c1 = c.AddElem(v1)
                        c2 = c1.AddElem(v2)
                        bisect.insort(returnlist, c2)
    return returnlist





#In [33]: c = Configuration([v0,v1,v2,v3,v4,v5])

#In [34]: c.Laplacian = [-1,1,-1,1,-1,1]

#In [35]: c.ObtainValue(10)
#Out[35]: [5.9776578812861292, 6.9878416372343634e-10]

###################################################
### Generation of lists of connected components ###
###################################################

# This list contains all connected components of size at most 6 with height bounded by 1 #

# cc_hex_list = RecursivelyBuildComponents()

# This list contains all connected components of size 7 with height bounded by 1
# hex_level_7 = GenNextLevel(cc_hex_list[5], 1)

# This list contains all connected components of size 8 with height bounded by 1
# hex_level_8 = GenNextLevel(hex_level_7[1], 1)


# hex_level_9 = GenNextLevel(hex_level_8[1], 1)

# In [15]: len(hex_level_9[1])
# Out[15]: 3

#In [17]: for c in hex_level_9[1]:
#    ...:     c.NumBins = 2
#    ...:     print c.ObtainValueCareful1()
#    ...:     
#7.79868782087
#8.50543874677
#8.14742703566

###############################################################################
####The above calc. shows that the largest configuration has size at most 8####
###############################################################################



################################################################################
# The best connected components with height at most 1 are enumerated as follows#
################################################################################

# In [352]: c4 = Obtain4C2Components(cc_list_hex[3])

#In [356]: for c in c4:
#     ...:     val = c.ObtainValueSimple(5)
#     ...:     if val < 6:
#     ...:         print [c, val]
#     ...: 

#In [354]: c6 = Obtain6C2Components(cc_list_hex[5])

#In [357]: for c in c6:
#     ...:     val = c.ObtainValueSimple(5)
#     ...:     if val < 6:
#     ...:         print [c, val]
#     ...:         
#[[[ 0.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 1.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]], 5.9763851811361848]


#In [18]: c8= Obtain8C2Components(hex_level_8[1])

#In [21]: for c in c8:
#    ...:     print c.ObtainValueSimple(2)
#    ...:     
#29.7570653118
#15.6295729165
#12.8762907506
#12.1812006772
#20.3124942868
#16.2451730212
#14.9277442017
#12.4547008605
#24.4176449444
#14.5859615349
#14.5859615349
#15.7802038445
#29.5583279304
#20.0315993274
#16.341713327
#16.2226042068
#16.9563465092
#13.6938211665
#17.2581377892
#13.403865247
#22.5397584277
#13.9036626982
#12.4430165124
#23.603331445
#9.909838888
#21.1365622573
#26.2560615943
#13.0725699876
#16.8344893203
#27.4813403682
#12.4244486799
#22.1286423405
#22.3193512369

##############################################################################

# cc2_list_hex = RecursivelyBuildComponents2()

# cc22_list_hex = RecursivelyBuildComponents22()
#In [350]: len(cc2_list_hex)
#Out[350]: 8

#In [351]: len(cc22_list_hex)
#Out[351]: 6


#################################################
#################################################
# This calculation checks the values in Lemma 50 #
#In [37]: v0 = Vertex(0,0,0)

#In [38]: v1 = Vertex(0,0,1)

#In [39]: v2 = Vertex(1,0,0)

#In [40]: c = Configuration([v0])

#In [41]: c.LaplaceConstraint = [2]

#In [42]: c.ObtainValueCareful2()
#/usr/lib/python2.7/dist-packages/scipy/optimize/_minimize.py:385: RuntimeWarning: Method SLSQP does not use Hessian information (hess).
#  RuntimeWarning)
#Out[42]: 3.4999999999978497

#In [43]: c1 = Configuration([v0, v1])

#In [44]: c1.LaplaceConstraint = [2,2]

#In [46]: c1.NumBins = 2

#In [47]: c1.ObtainValueCareful2()
#Out[47]: 3.9826034686186644

#In [48]: c2 = Configuration([v0,v2])

#In [49]: c2.LaplaceConstraint = [2,2]

#In [56]: c2.NumBins = 4

#In [57]: c2.ObtainValueCareful2()
#Out[57]: 5.9811984884674105

#In [60]: c.LaplaceConstraint = [1]

#In [62]: c.NumBins = 10

#In [63]: c.ObtainValueCareful1()
#Out[63]: 1.3537809676341894

##########

# This section checks the values in Lemma 52:

#In [58]: c1.LaplaceConstraint = [1,1]

#In [59]: c1.ObtainValueCareful1()
#Out[59]: 1.875000222691907

#In [64]: c2.LaplaceConstraint =  [1,1]

#In [69]: c2.NumBins = 10

#In [70]: c2.ObtainValueCareful1()
#Out[70]: 2.590553641113512

#In [71]: c1.LaplaceConstraint = [1,0]

#In [72]: c1.NumBins = 2

#In [73]: c1.ObtainSignedOptimizationCareful()
#Out[73]: 1.7205706431839667

#In [77]: l
#Out[77]: [(-1,0,1) , (0,-1,1) , (0,0,0) , (0,0,1) ]

#In [78]: c = Configuration(l)

#In [79]: c.LaplaceConstraint = [0,0,1,0]

#In [81]: c.NumBins = 2

#In [82]: c.ObtainSignedOptimizationCareful()
#Out[82]: 2.9296238684464337

#In [84]: l
#Out[84]: 
#[(-1,0,0) ,
# (-1,0,1) ,
# (-1,1,0) ,
# (0,-1,0) ,
# (0,-1,1) ,
# (0,0,0) ,
# (0,0,1) ,
# (0,1,0) ,
# (1,-1,0) ,
# (1,0,0) ]

#In [85]: c = Configuration(l)

#In [86]: c.LaplaceConstraint = [0,0,0,0,0,1,0,0,0,0]

#In [87]: c.NumBins = 1

#In [88]: c.ObtainSignedOptimizationCareful()
#Out[88]: 4.562746066994829





######################
# This section checks the values in Lemma 55 #

#In [28]: v0 = Vertex(0,0,0)

#In [29]: v1 = Vertex(0,0,1)

#In [30]: c = Configuration([v0,v1])

#In [31]: c.LaplaceConstraint = [1,1]

#In [32]: c.ObtainSignedOptimizationCareful()
#/usr/lib/python2.7/dist-packages/scipy/optimize/_minimize.py:385: RuntimeWarning: Method SLSQP does not use Hessian information (hess).
#  RuntimeWarning)
#Out[32]: 3.7955921367155465

# Q({0,v,v1}, \nu) \geq P({0,v1},1) \geq 2.59

########### This section reduces the number of connected components to be considered.  It is still necessary to eliminate the cases of 2's

#In [14]: for v in cc_hex_list[0]:
#    ...:     print v.value
#    ...:     
#    ...:     
#1.35424868903




# Thus there are at most 4 connected components.

#In [21]: for v1 in cc_hex_list[2]:
#    ...:     if v1.bestValue +3* 1.35 < target_value:
#    ...:         print v1
#    ...:         

#In [22]: for v1 in cc_hex_list[1]:
#    ...:     for v2 in cc_hex_list[1]:
#    ...:         if v1.bestValue + v2.bestValue + 2*1.35 < target_value:
#    ...:             print [v1, v2]
#    ...:             


# The last two calculations show that if there are 4 connected components, each is a singleton.

#########

# This verifies that there is not a component of size at least 6 which can be paired with two singletons

#In [24]: for v in cc_hex_list[5]:
#    ...:     if v.bestValue + 2*1.35 < target_value:
#    ...:         print v
#    ...:


# This enumerates the components of size 4 which can be paired with two singletons

#In [23]: for v in cc_hex_list[3]:
#    ...:     if v.bestValue + 2*1.35 < target_value:
#    ...:         print v
#    ...:         
#[[ 1.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]]
#[[ 1.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 1.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]]
#[[ 0.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 1.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]]

#######

# This calculation enumerates all configurations of size 2 and 3 which can be apired with a singleton

#In [25]: for v1 in cc_hex_list[1]:
#    ...:     for v2 in cc_hex_list[2]:
#    ...:         if v1.bestValue + v2.bestValue + 1.35 < target_value:
#    ...:             print [v1,v2]
#    ...:             
#[[[ 1.  0.  0.]
# [ 0.  1.  0.]
# [ 0.  0.  0.]], [[ 1.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]]]

# This calculation verifies that there is not a configuration of size 5 which can be paired with a configuration of size 2 and a singleton

#In [26]: for v1 in cc_hex_list[1]:
#    ...:     for v2 in cc_hex_list[4]:
#    ...:         if v1.bestValue + v2.bestValue + 1.35 < target_value:
#    ...:             print [v1,v2]
#    ...:

# This calculation verifies that there are not configurations of size 3 and 4 which can be paired with a singleton

#In [27]: for v1 in cc_hex_list[2]:
#    ...:     for v2 in cc_hex_list[3]:
#    ...:         if v1.bestValue + v2.bestValue + 1.35 < target_value:
#    ...:             print [v1,v2]
#    ...:             

# This calculation verifies that if there are three components of size 2, then each is a pair of adjacent vertices.

#In [28]: for v1 in cc_hex_list[1]:
#    ...:     for v2 in cc_hex_list[1]:
#    ...:         for v3 in cc_hex_list[1]:
#    ...:             if v1.bestValue + v2.bestValue + v3.bestValue < target_valu
#    ...: e:
#    ...:                 print [v1,v2,v3]
#    ...:                 
#    ...:             
#[[[ 1.  0.  0.]
# [ 0.  1.  0.]
# [ 0.  0.  0.]], [[ 1.  0.  0.]
# [ 0.  1.  0.]
# [ 0.  0.  0.]], [[ 1.  0.  0.]
# [ 0.  1.  0.]
# [ 0.  0.  0.]]]

## This calculation checks that if there are three pairs of adjacent vertices which are 3-connected, then the value is too large
 #(First generate all sets of two pairs at distance 3)
# In [300]: l = TwoPairsAdjacent()

# (Next add on a third pair at distance 3)

#In [307]: l1 = size4adjacentpairs(l)

# (Next assign values in C2)

# In [310]: l2 = Obtain6C2Components(l1)

# (Now calculate values of xi near 0 and estimate f(xi))

#In [312]: for c in l2:
#     ...:     print c.ObtainValueSimple(2)
#     ...:     
#10.6437184307
#10.4551640395
#12.814205331
#13.3812417216
#11.8048843461
#11.8048843461
#11.6392773525
#15.758727969
#15.758727969
#15.7640944976
#15.7640944976
#8.93608082106
#11.8048843461
#11.8048843461
#15.7640944976
#15.7640944976
#16.4953772694
#16.4953772694
#12.814205331
#10.6437184307
#13.3812417216
#11.8048843461
#11.8048843461
#11.6392773525
#10.4551640395
#11.8048843461
#11.8048843461
#8.93608082106



# This calculation rules out two size 2 configurations with a size 4 configuration

#In [29]: for v1 in cc_hex_list[1]:
#    ...:     for v2 in cc_hex_list[1]:
#    ...:         for v3 in cc_hex_list[3]:
#    ...:             if v1.bestValue + v2.bestValue + v3.bestValue < target_valu
#    ...: e:
#    ...:                 print [v1,v2,v3]
#    ...:                 
#    ...:     

# This calculation rules out two size 3 configurations with a size 2 configuration

#In [30]: for v1 in cc_hex_list[1]:
#    ...:     for v2 in cc_hex_list[2]:
#    ...:         for v3 in cc_hex_list[2]:
#    ...:             if v1.bestValue + v2.bestValue + v3.bestValue < target_valu
#    ...: e:
#    ...:                 print [v1,v2,v3]
#    ...:                 
#    ...:  




#######

# The value of Lemma 57 is checked here

#In [28]: l
#Out[28]: 
#[(0,0,0) ,
# (0,0,1) ,
# (0,1,0) ,
# (0,1,1) ,
# (1,0,0) ,
# (1,0,1) ,
# (1,1,0) ,
# (1,1,1) ]

#In [29]: c = Configuration(l)

#In [30]: c.LaplaceConstraint = [0,1,0,0,0,0,-1,0]

#In [31]: c.ObtainSignedOptimizationCareful()
#Out[31]: 5.033936715896039



###########
# This rules out pairing a component of size 7 with another component
#In [38]: hex_level_7 = GenNextLevel(cc_list_hex[5])

#In [47]: for v in hex_level_7[1]:
#    ...:     bv = min(bv, v.bestValue)
#    ...:     

#In [48]: bv
#Out[48]: 5.002386543209233




# There is only one configuration of size 6 which can be paired with a configuration of size 2, and the configuration of size 2 must be two adjacent nodes

#In [6]: for c in hex_list[5]:
#   ...:     if c.bestValue + 1.87 < 6:
#   ...:         print c
#   ...:         
#[[ 0.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 1.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]]

#In [7]: for c in hex_list[5]:
#   ...:     if c.bestValue + 2.48 < 6:
#   ...:         print c
#   ...:  

#This shows that there are no pairs of configurations with one of size 3 and the other of size 5 below the minimum

#In [14]: for c1 in hex_list[2]:
#    ...:     for c2 in hex_list[4]:
#    ...:         if c1.bestValue + c2.bestValue < 6:
#    ...:             print [c1, c2]
#    ...:  

####################


#In [63]: l = []

#In [64]: for c in hex_list[4]:
#    ...:     if c.bestValue + 1.35 < 6:
#    ...:         l1 = FindSize5SingletonPair(c)
#    ...:         for c1 in l1:
#    ...:             l.append(c1)
#    ...:
#In [66]: len(l)
#Out[66]: 36

# This computation proves that there is not a 5-1 pair which obtains the optimum

#In [68]: for c in l:
#    ...:     out = c.ObtainValue(5)
#    ...:     if out[0]-out[1] < 6:
#    ...:         print c
#    ...:         


####################


# The values in Lemma 58 are verified here: #


#In [169]: l
#Out[169]: [(-1,0,1) , (0,-1,1) , (0,0,0) , (0,0,1) , (0,1,0) , (1,0,0) ]

#In [170]: c1 = Configuration(l)
#In [177]: c1.NumBins = 2
#In [179]: c1.LaplaceConstraint = [0,0,1,-1,0,0]

#In [180]: c1.ObtainSignedOptimizationCareful()
#Out[180]: 2.99541719025144

#In [157]: l1
#Out[157]: [(-1,0,1) , (0,-1,1) , (0,0,0) , (0,0,1) , (1,-1,1) , (1,0,0) , (1,0,1) ]

#In [158]: c = Configuration(l1)

#In [159]: c.LaplaceConstraint = [0,0,1,0,0,-1,0]

#In [160]: c.ObtainSignedOptimizationCareful()
#Out[160]: 4.476995712507636


### v0 = (0,0,0), v2 = (1,0,0)
##In [161]: c = Configuration([v0, v2])

#In [162]: c.LaplaceConstraint = [1,-1]

#In [164]: c.NumBins = 2

#In [165]: c.ObtainSignedOptimizationCareful()
#Out[165]: 2.761292960040067

#In [157]: l1
#Out[157]: [(-1,0,1) , (0,-1,1) , (0,0,0) , (0,0,1) , (1,-1,1) , (1,0,0) , (1,0,1) ]
#In [205]: c.LaplaceConstraint = [0,0,1,0,0,1,0]

#In [206]: c.ObtainSignedOptimizationCareful()
#Out[206]: 7.209401545755421


#####################################

# This rules out 4-4 pairs #
# In [285]: bv = 10

# In [286]: for c in cc_list_hex[3]:
#     ...:     c.NumBins = 1
#     ...:     bv = min(bv, c.ObtainValueCareful1())
#     ...:     

# In [287]: bv
# Out[287]: 3.2687925964862368


########

# Next, 4-2 pairs are considered.  The calculations above show that if a 4-2 pair occurs, the distance between the pairs is exactly 3.

## This rules out a component of size 4 with a component of 2 adjacent nodes

#In [249]: for c in cc_list_hex[3]:                
#     ...:     if c.ObtainValueCareful1() + 1.87 < 6:
#     ...:         l.append(c)  

#In [240]: l0 = size4adjacentpairs(l)

#In [241]: len(l0)
#Out[241]: 242

#In [242]: l0_C2 = Obtain6C2Components(l0)

#In [243]: len(l0_C2)
#Out[243]: 34

#In [244]: l_ = []

#In [245]: for c in l0_C2:
#     ...:     if c.ObtainValueSimple(3) < 6:
#     ...:         l_.append(c)
#     ...:         

#In [246]: len(l_)
#Out[246]: 0


# This rules out a component of size 4 together with a component of size 2  with nodes at distance 2 from each other. #

# In [252]: l1 = []

# In [253]: for c in cc_list_hex[3]:
#     ...:     if c.ObtainValueCareful1() + 2.59 < 6:
#     ...:         bisect.insort(l1, c)
#     ...:    

# In [256]: l10 = size4dist2pairs(l1)

# In [258]: l10_ = Obtain6C2Components(l10)

# In [261]: for c in l10_:
#     ...:     print c.ObtainValueSimple(3)
#     ...:     
#27.3341950676
#25.0442963503
#28.0708536197
#21.0964525548
#26.5113658559
#21.27703421
#21.27703421
#22.4512407321
#26.9739176001
#22.0547150197
#22.0547150197
#29.5553887225
#26.6676610089
#28.485342848
#31.7972222113
#32.4706319332


###### The set of pairs of components of size 3 is narrowed as follows #

# In [271]: for c in cc_list_hex[2]:
#     ...:     c.NumBins = 4
#     ...:     print c.ObtainValueCareful1()
#     ...:     
# 2.62826448834
# 3.14211165157
# 3.82373274251
# 3.27604963368
# 3.8237327514
# 3.14211207589
# 3.70722137002

# In [272]: cc_list_hex[2][0].vertices
# Out[272]: [(0,0,0) , (0,0,1) , (0,1,0) ]

# In [273]: cc_list_hex[2][1].vertices
# Out[273]: [(0,0,0) , (0,0,1) , (0,1,1) ]

# In [274]: cc_list_hex[2][3].vertices
# Out[274]: [(0,0,0) , (0,1,0) , (1,0,0) ]

# In [275]: cc_list_hex[2][5].vertices
# Out[275]: [(0,0,1) , (0,1,0) , (1,0,1) ]

#In [278]: c = cc_list_hex[2][3]

#In [279]: c
#Out[279]: 
#[[ 1.  0.  0.  1.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 1.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]]

#In [280]: c.LaplaceConstraint = [1,1,-1]

#In [284]: c.ObtainSignedOptimizationCareful()
#Out[284]: 3.841789659525243


######

# The case of two pairs of components of size 2 is eliminated as follows:


#### The following calculation eliminates the case of two pairs at distance 2 from each other

#In [122]: l = findDist2Pairs()
#In [124]: l1 = Obtain4C2Components(l)

#In [131]: for c in l1:
#     ...:     print c.ObtainValueSimple(5)
#     ...:     

#31.4255485404
#28.2510995637
#26.5286666629
#19.746127746
#24.9699483447
#14.3200275786
#17.3494696187
#19.8333258682
#17.2719397717



#### The following calculation eliminates the case of two pairs of adjacents nodes

#In [31]: l = findAdjacentPairs()
#In [33]: l1 = Obtain4C2Components(l)
#In [37]: for c in l1:
#    ...:     print c.ObtainValueSimple(5)
#    ...:     
#    ...:     
#15.1666084083
#11.930071881
#11.2803024973
#14.8173375929
#11.9300718811
#15.1666084083
#15.1666084083
#11.9300718811
#14.8173375929
#11.2803024973
#11.930071881
#15.1666084083


########## This concludes the consideration of configurations having height at most 1


#@#@#@#@

# The values in Lemma 59 are verified here #

#In [336]: c1.vertices
#Out[336]: [(0,0,0) , (0,0,1) ]

#In [337]: c1.LaplaceConstraint = [2,2]

#In [338]: c1.ObtainSignedOptimizationCareful()
#Out[338]: 12.0

#In [339]: c1.LaplaceConstraint = [2,1]

#In [340]: c1.ObtainSignedOptimizationCareful()
#Out[340]: 7.817010499779817


#######


# Next consider the case of a single node of height 2.  There are at most 2 components in this case.  Consider first the case of two components.

# The following calculation shows that there cannot be a component of size at least 6 having a single 2, with a second component

#In [117]: for c in cc2_hex_list[5]:
#     ...:     if c.bestValue + 1.35 < 6:
#     ...:         print [c, c.bestValue]
#     ...:         
#     ...:    

# If there is a component with a single 2, of size 5, then it must be the following 1.  However, it would have to be paired with a component of size 2, and the P value is too large

#In [116]: for c in cc2_hex_list[4]:
#     ...:     if c.bestValue + 1.35 < 6:
#     ...:         print [c, c.bestValue]
#     ...:         
#     ...:     
#[[[ 0.  0.  0.  1.  0.  0.]
# [ 0.  0.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]], 4.543792605647258]

## Now consider components of size 4, having a single 2 node.  Those which may be paired with a singleton are as follows: #

#In [115]: for c in cc2_hex_list[3]:
#     ...:     if c.bestValue + 1.35 < 6:
#     ...:         print [c, c.bestValue]
#     ...:         
#     ...:     
#[[[ 0.  0.  0.  1.  0.  0.]
# [ 0.  0.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]], 4.50260817387589]
#[[[ 1.  0.  0.]
# [ 0.  1.  0.]
# [ 0.  0.  0.]
# [ 1.  0.  0.]
# [ 0.  1.  0.]
# [ 0.  0.  0.]], 4.502608125616599]
#[[[ 0.  0.  0.  1.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]], 4.543792649063275]
#[[[ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]], 3.5]


# There is a single component with a 2 on 3 nodes which can be paired with a pair of nodes which are adjacent to each other, as the following calculation shows:

#In [323]: for c in cc2_list_hex[2]:
#     ...:     if c.ObtainValueCareful2() < 6-1.87:
#     ...:         print c
#     ...:         
#[[ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]]


## The above configurations were analyzed, and this reduces to considering single connected components.

# There is only 1 component of size 3 which can be given a Laplacian assignment in C2

#In [328]: c.vertices
#Out[328]: [(-1,0,0) , (0,0,0) , (1,0,0) ]

#In [329]: c.Laplacian = [-1,2,-1]

#In [331]: c.ObtainValueSimple(3)
#Out[331]: 11.86590153891725

#In [332]: c.LaplaceConstraint = [-1,2,-1]

#In [333]: c.ObtainSignedOptimizationCareful()
#Out[333]: 6.393271261338805



# The size 5 components are as follows


#In [114]: l = Obtain5C2_2(cc2_hex_list[4])
#In [114]: len(l)
#Out[114]: 7

#### The last configuration is equivalent to the extremal `benzene ring' configuration

#15.1666084083
#30.2278269574
#12.1443460865
#28.7521947447
#14.8173375929
#34.1292238252
#5.97638517989

#In [110]: l[6]
#Out[110]: 
#[[ 0.  0.  0.  1.  0.  0.]
# [ 0.  1.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 1.  0.  0.  1.  0.  0.]
# [ 0.  0.  0.  0.  1.  0.]
# [ 0.  0.  0.  0.  0.  0.]]

#In [111]: l[6].vertices
#Out[111]: [(-1,-1,1) , (-1,0,0) , (0,-1,0) , (0,0,0) , (0,0,1) ]

#In [112]: l[6].ObtainValue(5)
#Out[112]: [5.9776577502107227, 1.3498156646490999e-07]

#In [113]: l[6].Laplacian
#Out[113]: array([ 1., -1., -1.,  2., -1.])

##### This component is equivalent to the extremal configuration with height one nodes around a hexagon #

# There were no size 7 components to consider

#In [116]: l3 = Obtain7C2_2(cc2_hex_list[6])

#In [117]: len(l3)
#Out[117]: 0


#@#@#@#@#@#@#@#@#@#@#@#@#@

# This has reduced to the case of two adjacent nodes of height 2, which necessarily have opposing sign.


# First consider the case of two connected components.  The following calculation rules out the case that the component with 2 2's has size 4 or more

#In [345]: for c in cc22_list_hex[2]:
#     ...:     print c.ObtainValueCareful2()
#     ...:     
#5.78403368534
#5.53892651767
#5.53892651771
#5.78403350411
#5.78403350769
#5.53892651763
#5.18824628196
#5.58090960918
#5.55209004077
#5.55209004077
#5.18824628232
#5.49777359996
#4.80549681652
#4.66265833703
#5.58090948986
#5.55209004743
#4.66265833703
#5.743481656

# If the component with 2 2's has size 3, it can only be of the following form:

#In [348]: for c in cc22_list_hex[1]:
#     ...:     print [c, c.ObtainValueCareful2()]
#     ...:     
#[[[ 1.  0.  0.]
# [ 0.  0.  0.]
# [ 0.  0.  0.]
# [ 1.  0.  0.]
# [ 0.  1.  0.]
# [ 0.  0.  0.]], 5.1802719953840874]
#[[[ 0.  0.  0.]
# [ 0.  1.  0.]
# [ 0.  0.  0.]
# [ 1.  0.  0.]
# [ 0.  1.  0.]
# [ 0.  0.  0.]], 4.328163528964978]
#[[[ 0.  0.  0.  1.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]
# [ 1.  0.  0.  0.  0.  0.]
# [ 0.  1.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  0.  0.]], 5.180271995382391]

#In [349]: cc22_list_hex[1][1].vertices
#Out[349]: [(-1,0,1) , (0,0,0) , (0,0,1) ]

# A connected component would have to be a singleton, which places it at location v2.


# The cases of a single connected component with 2 2's are ruled out as follows #

#In [134]: l = Obtain4C2_22(cc22_hex_list[2])

#In [135]: l
#Out[135]: 
#[[[ 0.  0.  0.  0.  0.  0.]
#  [ 0.  1.  0.  0.  0.  0.]
#  [ 0.  0.  0.  0.  0.  0.]
#  [ 1.  0.  0.  1.  0.  0.]
#  [ 0.  1.  0.  0.  0.  0.]
#  [ 0.  0.  0.  0.  0.  0.]]]

#In [136]: for c in l:
#     ...:     print c.ObtainValue(5)
#     ...:     
#[8.9904709769146098, 0.00039461093960291494]

#In [149]: Obtain6C2_22(cc22_hex_list[4])
#Out[149]: []

#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$








