###############################################################################################################
# This file contains the implementation of the hex tiling configuration class.  A hex tiling configuration    #
# consists of a list of vertices in the triangular lattice, together with a prevector which assigns integer   #
# values to the nodes.  The class includes as methods geometry functions (translate, rotate, obtain neighbors)#
# as well as implementations of the Fourier integral to obtain the corresponding function xi, and optimization#
# routines for P, Q, P_j, Q_j. Usage of this file appears in the comment section at the bottom, which verifies#
# the corresponding section of spectral_gap_comp.pdf in the repository.                                       #
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



#######################################################################
# The value of gamma_hex, with an error estimate, is calculated here. #
#######################################################################

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



##################################################################################
# The best connected components with height at most 1 are enumerated as follows. #  
# Such a connected component has an even number of nodes, and is in C2.          #
##################################################################################

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

####################################################################################################
# Next build a list of connected components that can have one, or two adjacent, nodes of height 2. #
####################################################################################################

# cc2_list_hex = RecursivelyBuildComponents2()

# cc22_list_hex = RecursivelyBuildComponents22()
#In [350]: len(cc2_list_hex)
#Out[350]: 8

#In [351]: len(cc22_list_hex)
#Out[351]: 6


##################################################
# This calculation checks the values in Lemma 21 #
##################################################

#In [37]: v0 = Vertex(0,0,0)

#In [38]: v1 = Vertex(0,0,1)

#In [39]: v2 = Vertex(1,0,0)

#In [40]: c = Configuration([v0])

#In [41]: c.LaplaceConstraint = [2]

#In [42]: c.ObtainValueCareful2()
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

###############################################
# This section checks the values in Lemma 23: #
###############################################

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





##############################################
# This section checks the values in Lemma 26 #
##############################################

#In [28]: v0 = Vertex(0,0,0)

#In [29]: v1 = Vertex(0,0,1)

#In [30]: c = Configuration([v0,v1])

#In [31]: c.LaplaceConstraint = [1,1]

#In [32]: c.ObtainSignedOptimizationCareful()
#Out[32]: 3.7955921367155465

# Q({0,v,v1}, \nu) \geq P({0,v1},1) \geq 2.59

####################################################################################################
# This section reduces the number of connected components to be considered.  It is still necessary #
# to eliminate the cases of nodes of height 2, which is considered below.                          #
####################################################################################################

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

#########################################
# The value of Lemma 27 is checked here #
#########################################

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

#############################################
# The values in Lemma 28 are verified here: #
#############################################

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

##############################################################################
# This concludes the consideration of configurations having height at most 1 #
##############################################################################

#@#@#@#@

############################################
# The values in Lemma 29 are verified here #
############################################

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









