###############################################################################################################
# This file contains the implementation of the hex tiling configuration class.  A hex tiling configuration    #
# consists of a list of vertices in the triangular lattice, together with a prevector which assigns integer   #
# values to the nodes.  The class includes as methods geometry functions (translate, rotate, obtain neighbors)#
# as well as implementations of the Fourier integral to obtain the corresponding function xi, and optimization#
# routines for P, Q, P_j, Q_j. Usage of this file is contained in connected_component_hex.py.                 #
###############################################################################################################


import numpy as np
from scipy import signal
from scipy import misc
import bisect
from copy import copy, deepcopy
from scipy.optimize import minimize
import scipy.integrate as integrate
from numpy import sin, cos, pi
from scipy.integrate import dblquad
import cmath
import matplotlib.pyplot as plt
import itertools as it
import math
from scipy.optimize import linprog

target_value = 5.9776

########################################################################################################
# The next few functions provide a method of integrating in polar coordinates to avoid the singularity #
# of the Green's function at 0. Since the Fourier integrals are over (R/Z)^2, a spherical neighborhood #
# of 0 is integrated in polar coordinates, then the remainder of the integral is performed over the    #
# torus.                                                                                               #
########################################################################################################

# Bump function which is 1 on [0,1/2] and 0 and [1, infty]

def BumpFunction(x):
    def cutoff(t):
        if t <= 0:
            return 0
        else:
            return np.exp(-1./t)
    def zero_one_transition(t):
        return cutoff(t)/(cutoff(t) + cutoff(1-t))
    return zero_one_transition(-2*x + 2)

# Integration in polar coordinates for integrand supported in r<=.5

def PolarCoordinateIntegrate(integrand):
    def polar_integrand(r,thet):
        if r == 0:
            return 0
        else:
            return integrand(r*np.cos(thet), r*np.sin(thet))*r
        
    return integrate.nquad(lambda r, thet: polar_integrand(r,thet), [[0, .5], [0, 2*np.pi]])[0]

# Integration over the torus (R/Z)^2 with singularity at 0

def SingularityIntegral(integrand):
    def smooth_integrand(x,y):
        if x**2 + y**2 == 0:
            return 0
        else:
            return integrand(x,y)* (1-BumpFunction(4*(x**2 + y**2)))
    def singular_integrand(x,y):
        return integrand(x,y) * BumpFunction(4*(x**2 + y**2))
    SmoothIntegral = integrate.nquad(lambda x,y: smooth_integrand(x,y), [[-.5,.5],[-.5,.5]])[0]
    SingularIntegral = PolarCoordinateIntegrate(singular_integrand)
    return SmoothIntegral + SingularIntegral



# Cosine difference in optimization functional

def CosDif(x):
    return 1 - np.cos(2 *np.pi * x)

# The objective function is linearized by splitting $[1/4, 1/2]$ into bins and interpolating CosDif linearly between the endpoints

def ObjFnSafe(x, bins):
    y = abs(x)
    if y < 0.25:
        return CosDif(y)
    elif y == 0.5:
        return 2
    else:
        j = math.floor(4*bins *(y-0.25))
        lower_pt = 0.25 + j*(0.25/bins)
        upper_pt = 0.25 + (j+1)*(0.25/bins)
        slope = (CosDif(upper_pt) - CosDif(lower_pt))*4*bins
        return CosDif(lower_pt) + slope *(y-lower_pt)

# The derivative of the linearized objective function

def DerivSafe(x, bins):
    y = abs(x)
    if y < 0.25:
        return 2*np.pi * np.sin(2*np.pi*x)
    else:
        if y == 0.5:
            j = bins-1
        else:
            j = math.floor(4*bins*(y-.25))
        lower_pt = 0.25 + j*0.25/bins
        upper_pt = 0.25 + (j+1)*0.25/bins
        slope = (CosDif(upper_pt) - CosDif(lower_pt))*4*bins
        if x>0:
            return slope
        else:
            return -slope

# The second derivative of the linearized objective function

def SecondDerivSafe(x):
    if (abs(x) < 0.25):
        return 4*np.pi**2 *np.cos(2*np.pi*x)
    else:
        return 0

#######################################
# Helper functions to work with lists #
#######################################

# Takes a list and returns it in sorted order

def SortList(l):
    returnlist = []
    for e in l:
        bisect.insort(returnlist,e)
    return returnlist

# Tests whether elem is present in a sorted list

def NotPresent(sort_list, elem):
    index = bisect.bisect(sort_list, elem)
    if index == 0 or elem != sort_list[index-1]:
        return True
    return False

# Two sorted lists are merged into a single sorted list

def MergeList(l1, l2):
    returnlist = SortList(l1)
    for e in l2:
        if NotPresent(returnlist, e):
	    bisect.insort(returnlist,e)
    return returnlist

###################################################################################################
# This class implements vertices in the triangular lattice. The coordinates give the location of  #
# x v_1 + y v_2 + parity v.  This object can be sorted in a sorted list, and can locate its       #
# neighbors in the lattice graph                                                                  #
###################################################################################################

class Vertex(object):
    def __init__(self, x_coor, y_coor, parity):
        self.x_coor = x_coor
        self.y_coor = y_coor
	self.parity = parity
    def __lt__(self, other):
        if self.x_coor < other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor < other.y_coor:
            return True 
	elif self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.parity < other.parity:
	    return True
        else:
            return False
    def __le__(self, other):
        if self.x_coor < other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor < other.y_coor:
            return True
	elif self.x_coor == other.x_coor and self.y_coor== other.y_coor and self.parity <= other.parity:
	    return True
        else:
            return False
    def __gt__(self, other):
        if self.x_coor > other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor > other.y_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.parity > other.parity:
            return True     
   	else:
            return False
    def __ge__(self, other):
        if self.x_coor > other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor > other.y_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.parity >= other.parity:
            return True
        else:
            return False
    def __eq__(self, other):
        if self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.parity == other.parity:
            return True
        else:
            return False
    def __ne__(self, other):
        if self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.parity == other.parity:
            return False
        else:
            return True
    def __getitem__(self, key):
        if key == 1:
            return self.x_coor
        if key == 2:
            return self.y_coor
	if key == 3:
	    return self.parity
        # Obtain neighbors of the vertex
    def ObtainNeighbors(self):
	if self.parity == 0:
	    v1 = Vertex(self.x_coor, self.y_coor, 1)
	    v2 = Vertex(self.x_coor-1, self.y_coor, 1)
	    v3 = Vertex(self.x_coor, self.y_coor-1, 1)
	    return SortList([v1, v2, v3])
	if self.parity == 1:
	    v1 = Vertex(self.x_coor, self.y_coor, 0)
	    v2 = Vertex(self.x_coor+1, self.y_coor, 0)
	    v3 = Vertex(self.x_coor, self.y_coor+1, 0)
	    return SortList([v1, v2, v3])
        # Obtains the distance 2 neighborhood of the vertex
    def ObtainDist2Neighbors(self):
        neighbors = self.ObtainNeighbors()
	returnlist = neighbors
	for n in neighbors:
	    returnlist = MergeList(returnlist, n.ObtainNeighbors())
	return returnlist
        # Obtains the distance 3 neighborhood of the vertex
    def ObtainDist3Neighbors(self):
        l = self.ObtainDist2Neighbors()
        l1 = l
        for v in l1:
            l = MergeList(l, v.ObtainNeighbors())
        return l
        # Obtains the neighbors at distance 1 in the triangular lattice of the vertex
    def ObtainTriangularLatticeNeighbors(self):
        v1 = Vertex(self.x_coor-1, self.y_coor, self.parity)
        v2 = Vertex(self.x_coor+1, self.y_coor, self.parity)
        v3 = Vertex(self.x_coor, self.y_coor-1, self.parity)
        v4 = Vertex(self.x_coor, self.y_coor+1, self.parity)
        v5 = Vertex(self.x_coor-1, self.y_coor+1, self.parity)
        v6 = Vertex(self.x_coor+1, self.y_coor-1, self.parity)
        return SortList([v1,v2,v3,v4,v5,v6])
        # Obtains the distance 2 neighborhood in the triangular lattice of the vertex
    def ObtainDist2TriangularLatticeNeighbors(self):
        l = self.ObtainTriangularLatticeNeighbors()
        returnlist = []
        for v in l:
            for w in v.ObtainTriangularLatticeNeighbors():
                if w != self and NotPresent(l, w) and NotPresent(returnlist, w):
                    bisect.insort(returnlist, w)
        return returnlist
    # Display
    def __str__(self):
        return "(" +str(self.x_coor) + "," +str(self.y_coor) + "," + str(self.parity) + ") "
    def __repr__(self):
        return "(" +str(self.x_coor) + "," +str(self.y_coor) + "," + str(self.parity) + ") "


        # Symmetries of the honeycomb tiling
    def Reflect(self):
        return Vertex(self.y_coor, self.x_coor, self.parity)

    def OrthogReflect(self):
        return Vertex(-self.y_coor, -self.x_coor, 1-self.parity)

    def CCR(self):
        return Vertex(-self.x_coor -self.y_coor - self.parity, self.x_coor, self.parity)

    def CCR2(self):
        return Vertex(self.y_coor, -self.x_coor - self.y_coor - self.parity, self.parity)

    def RCCR(self):
	return (self.CCR()).Reflect()

    def RCCR2(self):
	return (self.CCR2()).Reflect()

    def OR(self):
	return(self.Reflect()).OrthogReflect()

    def OCCR(self):
	return (self.CCR()).OrthogReflect()

    def OCCR2(self):
	return (self.CCR2()).OrthogReflect()

    def ORCCR(self):
	return (self.RCCR()).OrthogReflect()

    def ORCCR2(self):
	return (self.RCCR2()).OrthogReflect()

    def Translate(self, x, y):
        return  Vertex(self.x_coor+x, self.y_coor+y, self.parity)



####################################################################################################
# A configuration consists of a list of vertices in the support of a prevector.  The Laplacian is  #
# the function values at the points in the support.  This object can call an optimization program  #
# for either Q_j or P_j optimizations, or can use the Fourier integral to obtain the values of     #
# $\xi$ in a neighborhood of 0.                                                                    #
####################################################################################################

class Configuration(object):
    def __init__(self, vertices):
        self.error = -1
        self.vertices = SortList(vertices)
        self.n = len(vertices)
# This list contains constraint variables, and those variables belonging to two or more constraints
        self.large_variables = range(self.n)
# This is a list of nodes at distance 1 from the configuration
        self.neighbors = []
# This is a list of nodes at distance <= 2 from the configuration
        self.dist2neighbors = []
# The constraint for the P or Q optimization program
	self.LaplaceConstraint = np.ones(self.n)
# The values of the prevector on the configuration        
	self.Laplacian = np.zeros(self.n)
# Initialize neighbors
        for v in vertices:
            for n in v.ObtainNeighbors():
                if NotPresent(self.neighbors, n) and NotPresent(self.vertices, n):
                    bisect.insort(self.neighbors, n)
        for v in vertices:
            for n in v.ObtainDist2Neighbors():
                if NotPresent(self.dist2neighbors, n) and NotPresent(self.vertices, n):
                    bisect.insort(self.dist2neighbors, n)
        for i in range(len(self.neighbors)):
            v = self.neighbors[i]
            count = 0
            for v1 in v.ObtainNeighbors():
                if not NotPresent(self.vertices, v1):
                    count += 1
            if count >= 2:
                self.large_variables.append(self.n + i)
# Initialize bounds on coordinates
        self.x_max = vertices[0].x_coor
        self.x_min = vertices[0].x_coor
        self.y_max = vertices[0].y_coor
        self.y_min = vertices[0].y_coor
        for ve in vertices:
            self.x_max = max(self.x_max, ve.x_coor)
            self.x_min = min(self.x_min, ve.x_coor)
            self.y_max = max(self.y_max, ve.y_coor)
            self.y_min = min(self.y_min, ve.y_coor)
# Initialize an indicator function of the prevector, with value 1 at each vertex of the prevector
        self.prevector = np.zeros((3*(self.x_max - self.x_min + 1), 3*(self.y_max-self.y_min + 1)))
        for ve in vertices:
            self.prevector[3*(ve.x_coor-self.x_min)+ve.parity, 3*(ve.y_coor-self.y_min)+ve.parity] =1
        self.value = -1
        self.bestValue = -1
# The number of bins to use for optimizations P_j and Q_j
        self.NumBins = 1

# Display
    def __str__(self):
        return str(self.prevector)

    def __repr__(self):
        return str(self.prevector)

# Ordering
    def __lt__(self, other):
        if self.n < other.n:
            return True
        if other.n < self.n:
            return False
        iter1 = self.vertices.__iter__()
        iter2 = other.vertices.__iter__()
        while 1:
            try:
                f = iter1.next()
            except StopIteration:
                break
            s = iter2.next()
            if f < s:
                return True
            if s < f:
                return False
        return False
        
    def __le__(self, other):
        if self.n < other.n:
            return True
        if other.n < self.n:
            return False
        iter1 = self.vertices.__iter__()
        iter2 = other.vertices.__iter__()
        while 1:
            try:
                f = iter1.next()
            except StopIteration:
                break
            s = iter2.next()
            if f < s:
                return True
            if s < f:
                return False
        return True

    def __gt__(self, other):
        if self.n > other.n:
            return True
        if other.n > self.n:
            return False
        iter1 = self.vertices.__iter__()
        iter2 = other.vertices.__iter__()
        while 1:
            try:
                f = iter1.next()
            except StopIteration:
                break
            s = iter2.next()
            if f > s:
                return True
            if s > f:
                return False
        return False

    def __ge__(self, other):
        if self.n > other.n:
            return True
        if other.n > self.n:
            return False
        iter1 = self.vertices.__iter__()
        iter2 = other.vertices.__iter__()
        while 1:
            try:
                f = iter1.next()
            except StopIteration:
                break
            s = iter2.next()
            if f > s:
                return True
            if s > f:
                return False
        return True

    def __eq__(self, other):
        if (self.x_max - self.x_min == other.x_max - other.x_min) and (self.y_max - self.y_min == other.y_max - other.y_min):
            return (self.prevector == other.prevector).all()
        else:
            return False

    def __ne__(self, other):
        if (self.x_max - self.x_min == other.x_max - other.x_min) and (self.y_max - self.y_min == other.y_max - other.y_min):
            return (self.prevector != other.prevector).any()
        else:
            return True


    def setLaplaceConstraint(self, v):
	self.LaplaceConstraint = v

    def setLaplacian(self, v):
	self.Laplacian = v
        
# Returns a new configuration with the vertex v added to it
    def AddElem(self, v):
        newverts = []
        for v_ in self.vertices:
            v1 = Vertex(v_[1], v_[2], v_[3])
            newverts.append(v1)
        if NotPresent(newverts, v):
	    bisect.insort(newverts, v)
        return Configuration(newverts)

# This is used to build configurations by adding a vertex at distance 1 or 2 from the configuration
    def NextGeneration(self):
        gen = []
        for n in self.dist2neighbors:
            c = self.AddElem(n)
            cc = c.RotRefNormalize()
            if NotPresent(gen, cc):
                bisect.insort(gen, cc)
        return gen

# This is used to build configurations which have a single node of height 2    
    def NextGeneration2(self):
        v0 = Vertex(0,0,0)
        gen = []
        for n in self.dist2neighbors:
            c= self.AddElem(n)
            cc = c.RotRefQuotient()
            if NotPresent(gen, cc):
                i = bisect.bisect_left(cc.vertices, v0)
                cc.LaplaceConstraint[i] = 2
                bisect.insort(gen,cc)
        return gen

# This is used to build configurations which have a pair of adjacent nodes of height 2    
    def NextGeneration22(self):
        v0 = Vertex(0,0,0)
        v1 = Vertex(0,0,1)
        gen = []
        for n in self.dist2neighbors:
            c = self.AddElem(n)
            cc = c.ReflectionQuotient()
            if NotPresent(gen, cc):
                i = bisect.bisect_left(cc.vertices, v0)
                j = bisect.bisect_left(cc.vertices, v1)
                cc.LaplaceConstraint[i]=2
                cc.LaplaceConstraint[j]=2
                bisect.insort(gen,cc)
        return gen
    


# Some rigid motions of the plane to apply to configurations
    def Reflect(self):
        newverts = []
        for v in self.vertices:
            newverts.append(v.Reflect())
        return Configuration(SortList(newverts))

    def OrthogReflect(self):
	newverts = []
	for v in self.vertices:
	    newverts.append(v.OrthogReflect())
   	return Configuration(SortList(newverts))


    def CCR(self):
	newverts = []
	for v in self.vertices:
	    newverts.append(v.CCR())
	return Configuration(SortList(newverts))

    def CCR2(self):
	newverts = []
	for v in self.vertices:
 	    newverts.append(v.CCR2())
	return Configuration(SortList(newverts))

    def RCCR(self):
	return (self.CCR()).Reflect()

    def RCCR2(self):
	return (self.CCR2()).Reflect()

    def OR(self):
	return(self.Reflect()).OrthogReflect()

    def OCCR(self):
	return (self.CCR()).OrthogReflect()

    def OCCR2(self):
	return (self.CCR2()).OrthogReflect()

    def ORCCR(self):
	return (self.RCCR()).OrthogReflect()

    def ORCCR2(self):
	return (self.RCCR2()).OrthogReflect()

# Translate configurations by (x,y)

    def Translate(self, x, y):
        newverts = []
        for v in self.vertices:
            newverts.append(v.Translate(x,y))
        return Configuration(newverts)


# First vertex positive.  Configuration is shifted to have all x, y components non-negative, with some x and some y coordinate 0.
    def Normalize(self):
        cc = self.Translate(-self.x_min, -self.y_min)
        return cc

    def RotRefNormalize(self):
        cc1 = self.Normalize()
        cc2 = (self.CCR()).Normalize()
        cc3 = (self.CCR2()).Normalize()
        cc4 = (self.Reflect()).Normalize()
        cc5 = (self.RCCR()).Normalize()
        cc6 = (self.RCCR2()).Normalize()
        cc7 = (self.OrthogReflect()).Normalize()
        cc8 = (self.OR()).Normalize()
        cc9 = (self.OCCR()).Normalize()
        cc10 = (self.OCCR2()).Normalize()
        cc11 = (self.ORCCR()).Normalize()
        cc12 = (self.ORCCR2()).Normalize()
        l = SortList([cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8, cc9, cc10, cc11, cc12])
        return l[0]

# Quotient by reflections; these are symmetries which preserve a pair of nodes of height 2
    def ReflectionQuotient(self):
        cc1 = self
        cc2 = self.Reflect()
        cc3 = self.OrthogReflect()
        cc4 = self.OR()
        l = SortList([cc1,cc2,cc3,cc4])
        return l[0]

# Quotient by rotations and reflections which leave 0 fixed    
    def RotRefQuotient(self):
        cc1 = self
        cc2 = self.CCR()
        cc3 = self.CCR2()
        cc4 = self.Reflect()
        cc5 = self.RCCR()
        cc6 = self.RCCR2()
        l = SortList([cc1, cc2, cc3, cc4, cc5, cc6])
        return l[0]

# The number of variables in the P and Q optimization programs    
    def NumVariables(self):
        return len(self.vertices) + len(self.neighbors)

# The number of constraints in the P and Q optimization programs    
    def NumVertexVariables(self):
        return len(self.vertices)

# The objective function
    def ObtainFunctionValue(self, x):
        retval = self.NumVariables()
        for y in x:
            retval -= np.cos(2 * np.pi * y)
        return retval

# Linear interpolation of the objective function    
    def ObtainFunctionValueSafe(self,x):
        retval = 0
        for y in x:
            retval += ObjFnSafe(y, self.NumBins)
        return retval

# The derivative of the objective function
    def ObtainDerivative(self, x):
        deriv = np.zeros(self.NumVariables())
        for i in range(self.NumVariables()):
            deriv[i] = 2 * np.pi * np.sin(2 * np.pi *x[i])
        return deriv

# The derivative of the linearization of the objective function    
    def ObtainDerivativeSafe(self, x):
        deriv = np.zeros(self.NumVariables())
        for i in range(self.NumVariables()):
            deriv[i] = DerivSafe(x[i], self.NumBins)
        return deriv

# The Hessian of the objective function
    def ObtainHessian(self, x):
        Hessian = np.zeros((self.NumVariables(), self.NumVariables()))
        for i in range(self.NumVariables()):
            Hessian[i,i] = 4 * np.pi**2 * np.cos(2 * np.pi * x[i])
        return Hessian

# The Hessian of the linearization of the objective function
    def ObtainHessianSafe(self,x):
        Hessian = np.zeros((self.NumVariables(), self.NumVariables()))
        for i in range(self.NumVariables()):
            Hessian[i,i] = SecondDerivSafe(x[i], self.NumBins)
        return Hessian

# This version is used for the unsafe P programs
    def ObtainBounds(self):
        l=[]
        for i in range(self.NumVariables()):
            l.append((0, .5))
        return tuple(l)
        
# Bounds for the safe program P_j, bins = j
    def ObtainBoundsSafe(self, nl, vars, bins):
        lower_bound = [0,0,0,0,0]
        upper_bound = [0,0,0,0,0]
        for i in range(nl):
            lower_bound[i] = .25 + (.25 * bins[i])/self.NumBins
            upper_bound[i] = .25 + (.25 * (bins[i]+1.))/self.NumBins
        l = []
        for i in range(self.NumVariables()):
            if NotPresent(vars, i):
                l.append((0,.25))
            else:
                j = bisect.bisect(vars, i)-1
                l.append((lower_bound[j], upper_bound[j]))
        return tuple(l)

# Bounds for the unsafe program Q
    def ObtainSignedBounds(self):
        l = []
        for i in range(self.NumVariables()):
            l.append((-.5,.5))
        return tuple(l)

# Bounds for the safe program Q_j, j = bins
    def ObtainSignedBoundsSafe(self, nl, vars, bins, eps):
        upper_bound = [0,0,0,0,0]
        lower_bound = [0,0,0,0,0]
        for i in range(nl):
            lower_bound[i] = .25 + (.25*bins[i])/self.NumBins
            upper_bound[i] = .25 + (.25*(bins[i]+1.))/self.NumBins
            if eps[i] == -1:
                temp = - upper_bound[i]
                upper_bound[i]=-lower_bound[i]
                lower_bound[i] = temp
        l = []
        for i in range(self.NumVariables()):
            if NotPresent(vars, i):
                l.append((-.25, .25))
            else:
                j = bisect.bisect(vars, i)-1
                l.append((lower_bound[j], upper_bound[j]))
        return tuple(l)
    
    def ObtainConstraintMatrix(self):
        nv = self.NumVariables()
        nvv = self.NumVertexVariables()
        returnmat = np.zeros((nvv, nv))
        for i in range(nvv):
            returnmat[i,i] = 3
            for n in self.vertices[i].ObtainNeighbors():
                if not NotPresent(self.vertices, n):
                    ind = bisect.bisect(self.vertices, n)-1
                else:
                    ind = nvv + bisect.bisect(self.neighbors, n)-1
                returnmat[i, ind]=1
        return returnmat
            
# The unsafe optimization program P
    def ObtainOptimizationValue(self):
        x0 = np.zeros(self.NumVariables())
        for j in range(self.NumVariables()):
            x0[j] = .125        
        A = self.ObtainConstraintMatrix()
        def cons(x):
            return A.dot(x) - self.LaplaceConstraint
        constr = ({'type': 'ineq', 'fun': cons})
        res = minimize(lambda x: self.ObtainFunctionValue(x),x0, method ='SLSQP', bounds = self.ObtainBounds(), constraints =constr, jac = lambda x: self.ObtainDerivative(x), hess = lambda x: self.ObtainHessian(x) )
        if res.success:
            self.value = res.fun
            return   res.fun
        else:
            return -1
            
    def ObtainSignedConstraintMatrix(self):
        nv = self.NumVariables()
        nvv = self.NumVertexVariables()
        returnmat = np.zeros((nvv, nv))
        for i in range(nvv):
            returnmat[i,i] = 3
            for n in self.vertices[i].ObtainNeighbors():
                if not NotPresent(self.vertices, n):
                    ind = bisect.bisect(self.vertices, n)-1
                else:
                    ind = nvv + bisect.bisect(self.neighbors, n)-1
                returnmat[i, ind]=-1
        return returnmat

# The unsafe signed optimization program Q
    def ObtainSignedOptimizationValue(self):
        x0 = np.zeros(self.NumVariables())
        for j in range(self.NumVariables()):
            x0[j] = .125
        A = self.ObtainSignedConstraintMatrix()
        def cons(x):
            return A.dot(x) - self.LaplaceConstraint
        constr = ({'type': 'eq', 'fun': cons})
        res = minimize(lambda x: self.ObtainFunctionValue(x),x0, method ='SLSQP', bounds = self.ObtainSignedBounds(), constraints =constr, jac = lambda x: self.ObtainDerivative(x), hess = lambda x: self.ObtainHessian(x) )
        if res.success:
            self.value = res.fun
            return   res.fun
        else:
            return -1
        

# The safe optimization program P_j, to be used for prevectors with height <= 1.  This takes advantage of the speed-up that non-constraint variables which belong to only one constraint are bounded by 1/4.
    def ObtainValueCareful1(self):
        bestValue = self.ObtainOptimizationValue()
# Control the number of variables of size >= 1/4
        for num_large in range(5):
# Control which variable is large
            for vars in it.combinations(self.large_variables, num_large+1):
# Bin the range [1/4, 1/2]
                bin_iter = [it.product(range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins), range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins), range(self.NumBins), range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins), range(self.NumBins), range(self.NumBins), range(self.NumBins))]
                for bins in bin_iter[num_large]:
# This obtains a feasible initial point x0 
                    c = np.zeros(self.NumVariables())
                    c[0] = 1
                    b = self.LaplaceConstraint
                    A = self.ObtainConstraintMatrix()
                    opt_bounds = self.ObtainBoundsSafe(num_large+1, vars, bins)
                    res = linprog(c, A_eq = A, b_eq = b, bounds = opt_bounds)
                    if res.status == 0:                        
                        x0 = res.x
                        def cons(x):
                            return A.dot(x)-b
                        constr = ({'type': 'ineq', 'fun':cons})
# Perform the constrained minimization
                        res = minimize(lambda x: self.ObtainFunctionValueSafe(x), x0, method = 'SLSQP', bounds = opt_bounds, constraints = constr, jac = lambda x: self.ObtainDerivativeSafe(x), hess = lambda x: self.ObtainHessianSafe(x))
                        if res.success:
# Compare the answer to the answer from other bins
                            bestValue = min(bestValue, res.fun)
        self.bestValue = bestValue
        return self.bestValue

# The safe optimization program P_j, to be used for prevectors with a node of height 2.
    def ObtainValueCareful2(self):
        bestValue = self.ObtainOptimizationValue()
# Control the number of variables of size >= 1/4
        for num_large in range(5):
# Control which variables are large
            for vars in it.combinations(range(self.NumVariables()), num_large+1):
# Bin the large variables in subintervals of [1/4,1/2]
                bin_iter = [it.product(range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins), range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins), range(self.NumBins), range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins), range(self.NumBins), range(self.NumBins), range(self.NumBins))]
                for bins in bin_iter[num_large]:
# Obtain a feasible point x0 as initial guess
                    c = np.zeros(self.NumVariables())
                    c[0] = 1
                    b = self.LaplaceConstraint
                    A = self.ObtainConstraintMatrix()
                    opt_bounds = self.ObtainBoundsSafe(num_large+1, vars, bins)
                    res = linprog(c, A_eq = A, b_eq = b, bounds = opt_bounds)
                    if res.status == 0:
                        x0 = res.x
                        def cons(x):
                            return A.dot(x)-b
                        constr = ({'type': 'ineq', 'fun':cons})
# Perform the optimization program
                        res = minimize(lambda x: self.ObtainFunctionValueSafe(x), x0, method = 'SLSQP', bounds = opt_bounds, constraints = constr, jac = lambda x: self.ObtainDerivativeSafe(x), hess = lambda x: self.ObtainHessianSafe(x))
                        if res.success:
# Compare the optimization value to previous values
                            bestValue = min(bestValue, res.fun)
        self.bestValue = bestValue
        return self.bestValue

# This is the safe optimization program Q_j
    def ObtainSignedOptimizationCareful(self):
        bestValue = self.ObtainSignedOptimizationValue()
        lower_bound = [0,0,0,0,0]
        upper_bound = [0,0,0,0,0]
# Control the number of large variables
        for num_large in range(5):
# Control the sign of each large variable
            eps_ = [-1,1]
            eps_iter = [it.product(eps_), it.product(eps_,eps_), it.product(eps_,eps_,eps_), it.product(eps_,eps_,eps_,eps_), it.product(eps_,eps_,eps_,eps_,eps_)]
            for eps in eps_iter[num_large]:
# Bin the large variables size in a bin of [1/4, 1/2]
                for vars in it.combinations(range(self.NumVariables()), num_large+1):
                    bin_iter = [it.product(range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins), range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins), range(self.NumBins), range(self.NumBins)), it.product(range(self.NumBins), range(self.NumBins), range(self.NumBins), range(self.NumBins), range(self.NumBins))]  
                    for bins in bin_iter[num_large]:
# Check if the choice of bins automatically guarantees that the optimization value is too large
                        check = 0
                        for i in range(num_large+1):
                            lower_bound[i] = .25 + (.25 * bins[i])/self.NumBins
                            check += CosDif(lower_bound[i])
                        if check < target_value:
# Obtain a feasible initial guess x0
                            opt_bounds = self.ObtainSignedBoundsSafe(num_large+1, vars, bins, eps)
                            b = self.LaplaceConstraint
                            A = self.ObtainSignedConstraintMatrix()
                            c = np.zeros(self.NumVariables())
                            c[0]=1
                            res = linprog(c, A_eq = A, b_eq = b, bounds =opt_bounds)
                            if res.status ==0:
# Perform the minimization
                                x0 = res.x
                                def cons(x):
                                    return A.dot(x)-b
                                constr = ({'type': 'eq', 'fun':cons})
                                res = minimize(lambda x: self.ObtainFunctionValueSafe(x), x0, method = 'SLSQP', bounds = opt_bounds, constraints = constr, jac = lambda x: self.ObtainDerivativeSafe(x), hess = lambda x: self.ObtainHessianSafe(x))
                                if res.success:
# Compare the value to previous bin values
                                    val = res.fun
                                    bestValue = min(bestValue, val)
        self.bestValue = bestValue
        return self.bestValue

# Given the configuration, expands the configuration by adding all nodes at distance 1, with constraint value set to 0.  All 1/-1 assignments of the previous nodes are tested.    
    def ObtainOneBoundaryExpansionValue(self):
        l = MergeList(self.vertices, self.neighbors)
        c = Configuration(l)
        c.Laplacian = np.zeros(len(l))
        c.LaplaceConstraint = np.zeros(len(l))
        signs = it.product([-1,1], repeat = self.n-1)
        returnlist = []
        for eps in signs:
            count = 0
            for v in self.vertices:
                i = bisect.bisect(l, v)-1
                if count < self.n-1:
                    c.Laplacian[i] = eps[count]
                    c.LaplaceConstraint[i] = eps[count]
                else:
                    c.Laplacian[i] = 1
                    c.LaplaceConstraint[i] = 1
                count += 1
            returnlist.append([c, c.LaplaceConstraint, c.ObtainSignedOptimizationCareful()])
        return returnlist

# Given the configuration, expands the configuration by adding all nodes at distance at most 2, with constraint value set to 0. All 1/-1 assignments of the previous nodes are tested.
    def ObtainTwoBoundaryExpansionValue(self):
        l = MergeList(self.vertices, self.dist2neighbors)
        c = Configuration(l)
        c.Laplacian = np.zeros(len(l))
        c.LaplaceConstraint = np.zeros(len(l))
        signs = it.product([-1,1], repeat = self.n-1)
        returnlist = []
        for eps in signs:
            count = 0
            for v in self.vertices:
                i = bisect.bisect(l, v)-1
                if count < self.n-1:
                    c.Laplacian[i] = eps[count]
                    c.LaplaceConstraint[i] = eps[count]
                else:
                    c.Laplacian[i] = 1
                    c.LaplaceConstraint[i] = 1
                count += 1
            returnlist.append([c, c.LaplaceConstraint, c.ObtainSignedOptimizationCareful()])
        return returnlist

#This method obtains the Fourier integral representation of Xi in a neighborhood of 0, with n_1 v_1 + n_2 v_2 + n_3 v, |n_1|, |n_2| <= n
    
    def ObtainXiRepn(self, n):
        time_domain = np.zeros((2*n+1, 2*n+1,2))
        for a in range(2*n+1):
            for b in range(2*n+1):
                for p in range(2):
                    v = Vertex(a-n, b-n, p)
                    def integrand(x,y):
                        rvalue = 0
                        for k in range(self.n):
                            v_ = self.vertices[k]
                            if v_[3] == p:
                                rvalue += self.Laplacian[k] * np.cos(2 *np.pi * ((v_[1]-v[1])*x + (v_[2]-v[2])*y) )
                            else:
                                for v1 in v_.ObtainNeighbors():
                                    rvalue += (1./3.) *self.Laplacian[k] * np.cos(2*np.pi *((v1[1]-v[1])*x + (v1[2]-v[2])*y))
                        return rvalue / (2. - (2./3.) * (np.cos(2*np.pi * x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (x-y))))
                    time_domain[a,b,p] =  dblquad(lambda t,x: integrand(t,x), 0, 1, lambda x:0, lambda x:1, epsabs = 1.0e-11)[0]
        return time_domain

# This method calculates the L2 norm of Xi by Parseval
    
    def ObtainL2(self):
        def integrand1(x,y):
            rvalue1 = 0
            for k in range(self.n):
                v_ = self.vertices[k]
                if v_[3] == 0:
                    rvalue1 += self.Laplacian[k] * np.exp(2 *np.pi * 1j * (v_[1] * x + v_[2]*y))
                else:
                    for v1 in v_.ObtainNeighbors():
                        rvalue1 += (1./3.) *self.Laplacian[k] * np.exp(2*np.pi *1j*(v1[1]*x + v1[2]*y))
            return abs(rvalue1)**2/(2.-(2./3.) *(np.cos(2*np.pi * x) + np.cos(2*np.pi *y) + np.cos(2*np.pi*(x-y)) ))**2
        def integrand2(x,y):
            rvalue2 = 0
            for k in range(self.n):
                v_ = self.vertices[k]
                if v_[3] == 1:
                    rvalue2+= self.Laplacian[k] * np.exp(2*np.pi * 1j * (v_[1] * x + v_[2]*y))
                else:
                    for v1 in v_.ObtainNeighbors():
                        rvalue2 += (1./3.) * self.Laplacian[k] *np.exp(2*np.pi *1j * (v1[1]*x + v1[2]*y))
            return abs(rvalue2)**2/(2.-(2./3.)*(np.cos(2*np.pi*x) + np.cos(2*np.pi*y) + np.cos(2*np.pi*(x-y))))**2
        def integrand(x,y):
            return integrand1(x,y) + integrand2(x,y)
        return dblquad(lambda t,x: integrand(t,x),0,1,lambda x:0, lambda x:1, epsabs = 1.0e-11)

# This method calculates f(xi) by taking the values of xi in a neighborhood of 0 and combining with the L2 correction outside the neighborhood

    def ObtainValue(self, n):
        self.value = 0
        Xi = self.ObtainXiRepn(n)
        twonorm = self.ObtainL2()[0]
        approx2norm = 0
        for a in range(2*n+1):
            for b in range(2*n+1):
                for p in range(2):
                    self.value += 1 - np.cos(2 * np.pi * Xi[a,b,p])
                    approx2norm += Xi[a,b,p]**2
        remaining2norm = twonorm - approx2norm
        self.value += 2 * np.pi**2 * remaining2norm - (1./3.)*np.pi**4 * remaining2norm**2
        self.error = (1./3.)*np.pi**4 * remaining2norm**2
        return [self.value, self.error]

# This method calculates f(xi) without the L2 correction.
    
    def ObtainValueSimple(self, n):
        self.value = 0
        Xi = self.ObtainXiRepn(n)
        for a in range(2*n+1):
            for b in range(2*n+1):
                for p in range(2):
                    self.value += 1 - np.cos(2 * np.pi * Xi[a,b,p])
        return self.value


