############################################################################################################
# This file builds on configuration_hex.py to implement functions which recursively build connected        #
# component configurations, assign them C^2 values as prevectors, and combine several connected components.#
############################################################################################################

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
    counter = 0
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

# Recursively build all admissible components with nodes of height 1, currently this is stopped at configurations of size 6

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

# Given a list of configurations of size 6, tests all assignments of signs such that the resulting configuration is C2.  These calculations use that v = 1/3 (v_1 + v_2).

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

# Given a vector vx*v1 + vy*v2 + vz*v, returns the point in the hex tiling in standard form as a Vertex
def ObtainNormalizedVertex(vx, vy, vz):
    if vz %3 == 0:
        vx_= vx+ vz//3
        vy_= vy + vz//3
        vz_ = 0        
    elif vz%3 == 1:
        vx_ = vx + (vz-1)//3
        vy_ = vy + (vz-1)//3
        vz_ = 1
    v = Vertex(int(vx_), int(vy_), int(vz_))
    return v



# Given a configuration of size 3, iterates over all assignments of signs, and determines the location of a singleton such that the configuration is in C2



def FindSize3SingletonPair(c):
    returnlist = []
    for i_ in range(2):
        c1 = Configuration(c.vertices)
        c1.Laplacian = np.ones(3)
        c1.Laplacian[i_] = -1
        v1 = c1.vertices[i_]
	v = ObtainSingletonVertex(c1.vertices, c1.Laplacian)
        vx = 0
        vy = 0
        vz = 0
        for m in range(3):
            vx += c1.Laplacian[m] * c1.vertices[m].x_coor
            vy += c1.Laplacian[m] * c1.vertices[m].y_coor
            vz += c1.Laplacian[m] * c1.vertices[m].parity
        if vz % 3 != 2:
            v = ObtainNormalizedVertex(vx, vy, vz) 
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
            v = ObtainNormalizedVertex(vx, vy, vz) 
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



