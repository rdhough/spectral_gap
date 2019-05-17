############################################################################################################
# This file builds on configuration_triangular.py to implement functions which recursively build connected #
# component configurations, assign them C^2 values as prevectors, and combine several connected components.#
############################################################################################################

import numpy as np
from scipy import signal
from scipy import misc
import bisect
from copy import copy, deepcopy

target_value = 1.6936
target_value_1 = 1.6936-.44256
target_value_2 = 1.6936 - 0.6729

def GetCandidates(l):
    Candidates = []
    for c in l:
        for n in c.NextGeneration():
            if NotPresent(Candidates, n):
                bisect.insort(Candidates, n)
    return Candidates

# Starting from a configuration with a single vertex, configurations are built which can be obtained by adding a neighboring vertex.  If at some point the optimization value of a configuration is larger than the target value, then no futher configurations are built on top of the current one.

def RecursivelyBuildComponents():
    v = Vertex(0,0)
    c = Configuration([v])
    c.ObtainOptimizationValue()
    ConfigurationList = []
    ConfigurationList.append([c])
    while len(ConfigurationList[-1]) > 0:
        candidates = GetCandidates(ConfigurationList[-1])
        newlist = []
        for c in candidates:
            ov = c.ObtainValueCareful()
            if ov < target_value_2:
                bisect.insort(newlist, c)
        ConfigurationList.append(newlist)
    return ConfigurationList


#Given a list of configurations on 4 vertices, this finds ways of assigning +/- to the vertices to make them have moment 0 (this is the condition of being C^2).  

def Obtain4C2Components(l):
    returnlist = []
    for c in l:
        vx = 0
        vy = 0
        for v in c.vertices:
            vx += v[1]*.5
            vy += v[2]*.5
        v_ = Vertex(vx - c.vertices[0][1], vy-c.vertices[0][2])
        if not NotPresent(c.vertices, v_):
            ind = bisect.bisect(c.vertices, v_) - 1
            laplace = [-1,1,1,1]
            laplace[ind] = -1
            cc = Configuration(c.vertices)
            cc.setLaplacian(laplace)
            cc.setLaplaceConstraint(laplace)
            returnlist.append(cc)
    return returnlist

#Similar to above, but handles configurations with 6 vertices

def Obtain6C2Components(l):
    returnlist = []
    for c in l:
        vx = 0
        vy = 0
        for v in c.vertices:
            vx += v[1]*.5
            vy += v[2]*.5
        verts = c.vertices
        for i in range(1,5):
            for j in range(i+1,6):
                if (verts[0][1] + verts[i][1] + verts[j][1] == vx) and (verts[0][2] + verts[i][2] + verts[j][2] == vy):
                    laplace=[-1,1,1,1,1,1]
                    laplace[i] = -1
                    laplace[j] = -1
                    cc = Configuration(c.vertices)
                    cc.setLaplacian(laplace)
                    cc.setLaplaceConstraint(laplace)
                    returnlist.append(cc)
    return returnlist

# This function produces all pairs of adjacent vertices who have mutual distance 3.

def ObtainDist3PairNeighbors():
    v0 = Vertex(0,0)
    v1 = Vertex(1,0)
    l = []
    for v in v0.ObtainDist3Neighbors():
        v2 = Vertex(v.x_coor-1, v.y_coor)
        v3 = Vertex(v.x_coor+1, v.y_coor)
        if NotPresent( v0.ObtainDist2Neighbors(), v) and NotPresent( v1.ObtainDist2Neighbors(), v):
            if NotPresent( v0.ObtainDist2Neighbors(), v2) and NotPresent( v1.ObtainDist2Neighbors(), v2):
                c = Configuration([v0,v1,v,v2])
                if NotPresent(l, c):
                    bisect.insort(l,c)
            if NotPresent( v0.ObtainDist2Neighbors(), v3) and NotPresent( v1.ObtainDist2Neighbors(), v3):
                c1 = Configuration([v0,v1,v,v3])
                if NotPresent(l, c1):
                    bisect.insort(l, c1)
    return l

# This function checks for assignments of signs to a configuration of size 5 so that a singleton can be added at distance at least 3 and produce a configuration in C^2(T)

def ObtainSingletonC2_5(c):
    l = []
    for i in range(4):
        for offset in range(4-i):
            j = i + offset
            c.Laplacian = np.ones(5)
            c.Laplacian[i] = -1
            c.Laplacian[j] = -1
            x_0 = 0
            y_0 = 0
            for k in range(5):
                x_0 += c.Laplacian[k] * c.vertices[k].x_coor
                y_0 += c.Laplacian[k] * c.vertices[k].y_coor
            v = Vertex(x_0, y_0)
            v1 = Vertex(-x_0, -y_0)
            if NotPresent(c.dist2neighbors, v) and NotPresent(c.vertices, v):
                c1 = c.AddElem(v)
                l.append(c1)
            if NotPresent(c.dist2neighbors, v1) and NotPresent(c.vertices, v1):
                c2 = c.AddElem(v)
                l.append(c2)
    return l



