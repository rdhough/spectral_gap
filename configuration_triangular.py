###############################################################################################################
# This file contains the implementation of the triangular lattice configuration class.  A triangular          #
# lattice configuration consists of a list of vertices in the triangular lattice, together with a prevector   #
# which assigns integer values to the nodes.  The class includes as methods geometry functions (translate,    #
# rotate, obtain neighbors) as well as implementations of the Fourier integral to obtain the corresponding    #
# function xi, and optimization routines for P, Q, P_j, Q_j. Usage of this file is given in the comment       # 
# section at the bottom of the file, verifying the claims made in spectral_gap_comp.pdf from the repository.  #
###############################################################################################################

#import standard libraries

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

# Tests whether elem is present in a sorted list

def MergeList(l1, l2):
    returnlist = SortList(l1)
    for e in l2:
        if NotPresent(returnlist, e):
	    bisect.insort(returnlist,e)
    return returnlist

###################################################################################################
# This class implements vertices in the triangular lattice. The coordinates give the location of  #
# x v_1 + y v_2.  This object can be sorted in a sorted list, and can locate its neighbors in     #
# the lattice graph.                                                                              #
###################################################################################################
 
class Vertex(object):
    def __init__(self, x_coor, y_coor):
        self.x_coor = x_coor
        self.y_coor = y_coor
    def __lt__(self, other):
        if self.x_coor < other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor < other.y_coor:
            return True 
        else:
            return False
    def __le__(self, other):
        if self.x_coor < other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor <= other.y_coor:
            return True
        else:
            return False
    def __gt__(self, other):
        if self.x_coor > other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor > other.y_coor:
            return True
        else:
            return False
    def __ge__(self, other):
        if self.x_coor > other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor >= other.y_coor:
            return True
        else:
            return False
    def __eq__(self, other):
        if self.x_coor == other.x_coor and self.y_coor == other.y_coor:
            return True
        else:
            return False
    def __ne__(self, other):
        if self.x_coor == other.x_coor and self.y_coor == other.y_coor:
            return False
        else:
            return True
    def __getitem__(self, key):
        if key == 1:
            return self.x_coor
        if key == 2:
            return self.y_coor
        # Obtain neighbors of the vertex 
    def ObtainNeighbors(self):
        v1 = Vertex(self.x_coor-1, self.y_coor)
        v2 = Vertex(self.x_coor, self.y_coor-1)
        v3 = Vertex(self.x_coor, self.y_coor + 1)
        v4 = Vertex(self.x_coor+1, self.y_coor)
	v5 = Vertex(self.x_coor+1, self.y_coor-1)
	v6 = Vertex(self.x_coor-1, self.y_coor+1)
        return SortList([v1,v2,v3,v4,v5,v6])
        # Obtains the distance 2 neighborhood of the vertex 
    def ObtainDist2Neighbors(self):
	l = self.ObtainNeighbors()
        returnlist = l
	for e in l:
	    returnlist = MergeList(returnlist, e.ObtainNeighbors())
        return returnlist
        # Obtains the distance 3 neighborhood of the vertex 
    def ObtainDist3Neighbors(self):
        l = self.ObtainDist2Neighbors()
        l1 = l
        for v in l:
            l1 = MergeList(l1, v.ObtainNeighbors())
        return l1
        # Obtains the distance 4 neighborhood of the vertex
    def ObtainDist4Neighbors(self):
        l = self.ObtainDist2Neighbors()
        l1 = l
        for v in l:
            l1 = MergeList(l1, v.ObtainDist2Neighbors())
        return l1
	# Display
    def __str__(self):
        return "(" +str(self.x_coor) + "," +str(self.y_coor) + ") "
    def __repr__(self):
        return "(" +str(self.x_coor) + "," +str(self.y_coor) + ") "

####################################################################################################
# A configuration consists of a list of vertices in the support of a prevector.  The Laplacian is  #
# the function values at the points in the support.  This object can call an optimization program  #
# for either Q_j or P_j optimizations, or can use the Fourier integral to obtain the values of     #
# $\xi$ in a neighborhood of 0.                                                                    #
####################################################################################################
    
class Configuration(object):
    def __init__(self, vertices):
        self.vertices = SortList(vertices)
        self.n = len(vertices)
	self.LaplaceConstraint = np.ones(self.n)
	self.Laplacian = np.zeros(self.n)
        self.neighbors = []
        self.dist2neighbors = []
        self.dist3neighbors = []
        self.dist4neighbors = []
# Initialize a list of neighbors
        for v in vertices:
            for n in v.ObtainNeighbors():
                if NotPresent(self.neighbors, n) and NotPresent(self.vertices, n):
                    bisect.insort(self.neighbors, n)
        for v in vertices:
            for n in v.ObtainDist2Neighbors():
                if NotPresent(self.dist2neighbors, n) and NotPresent(self.vertices, n):
                    bisect.insort(self.dist2neighbors, n)
        for v in vertices:
            self.dist3neighbors = MergeList(self.dist3neighbors, v.ObtainDist3Neighbors())
            self.dist4neighbors = MergeList(self.dist4neighbors, v.ObtainDist4Neighbors())
# Initialize bounds on the coordinates
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
        self.prevector = np.zeros((self.x_max - self.x_min + 1, self.y_max-self.y_min + 1))
        for ve in vertices:
            self.prevector[ve.x_coor-self.x_min, ve.y_coor-self.y_min] =1
        self.value = -1
        self.bestValue = -1
# The number of bins to use for optimizations P_j and Q_j
        self.NumBins = 4

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

# Given a configuration, add a node at vertex v

    def AddElem(self, v):
        newverts = []
        for v_ in self.vertices:
            v1 = Vertex(v_[1], v_[2])
            newverts.append(v1)
        if NotPresent(newverts, v):
	    bisect.insort(newverts, v)
        return Configuration(newverts)

# Obtain all those configurations which may be obtained by adding a single connected vertex

    def NextGeneration(self):
        gen = []
        for n in self.dist2neighbors:
            c = self.AddElem(n)
            cc = c.RotRefNormalize()
            if NotPresent(gen, cc):
                bisect.insort(gen, cc)
        return gen


# Some rigid motions of the plane to apply to configurations
        

# Translate configurations by (x,y)

    def Translate(self, x, y):
        newverts = []
        for v in self.vertices:
            v_ = Vertex(v[1]+x, v[2]+y)
            newverts.append(v_)
        return Configuration(SortList(newverts))

# Rotations in 60 degree counterclockwise increments

    def CCR(self):
        newverts = []
        for v in self.vertices:
	    x_cor = - v[2]
            y_cor = v[1] + v[2]
            v1 = Vertex(x_cor, y_cor)
            newverts.append(v1)
        return Configuration(SortList(newverts))


    def CCR2(self):
        newverts = []
        for v in self.vertices:
	    x_cor = - (v[1]+v[2])
            y_cor = v[1] 
            v1 = Vertex(x_cor, y_cor)
            newverts.append(v1)
        return Configuration(SortList(newverts))


    def CCR3(self):
        newverts = []
        for v in self.vertices:
	    x_cor = - v[1]
            y_cor = - v[2]
            v1 = Vertex(x_cor, y_cor)
            newverts.append(v1)
        return Configuration(SortList(newverts))


    def CCR4(self):
        newverts = []
        for v in self.vertices:
	    x_cor =  v[2]
            y_cor = -(v[1] + v[2])
            v1 = Vertex(x_cor, y_cor)
            newverts.append(v1)
        return Configuration(SortList(newverts))


    def CCR5(self):
        newverts = []
        for v in self.vertices:
	    x_cor = v[1]+ v[2]
            y_cor = -v[1]
            v1 = Vertex(x_cor, y_cor)
            newverts.append(v1)
        return Configuration(SortList(newverts))

# Reflections in the dihedral group

    def Ref(self):
	newverts = []
	for v in self.vertices:
	    x_cor = v[1] + v[2]
	    y_cor = -v[2]
	    v1 = Vertex(x_cor, y_cor)
	    newverts.append(v1)
	return Configuration(SortList(newverts))

    def RCCR(self):
	c = self.Ref()
	return c.CCR()

    def RCCR2(self):
	c = self.Ref()
	return c.CCR2()

    def RCCR3(self):
	c = self.Ref()
	return c.CCR3()

    def RCCR4(self):
	c = self.Ref()
	return c.CCR4()

    def RCCR5(self):
	c = self.Ref()
	return c.CCR5()

    def Normalize(self):
        cc = self.Translate(-self.x_min, -self.y_min)
        return cc

# Quotient by the symmetry group

    def RotRefNormalize(self):
        cc = self.Normalize()
        cc1 = (self.CCR()).Normalize()
        cc2 = (self.CCR2()).Normalize()
        cc3 = (self.CCR3()).Normalize()
        cc4 = (self.CCR4()).Normalize()
        cc5 = (self.CCR5()).Normalize()
	rcc = (self.Ref()).Normalize()
	rcc1 = (self.RCCR()).Normalize()
	rcc2 = (self.RCCR2()).Normalize()
	rcc3 = (self.RCCR3()).Normalize()
	rcc4 = (self.RCCR4()).Normalize()
	rcc5 = (self.RCCR5()).Normalize()
        l = SortList([cc, cc1, cc2, cc3, cc4, cc5, rcc, rcc1, rcc2, rcc3, rcc4, rcc5])
        return l[0]

#This is the number of vertices in the support of the prevector, plus the number of vertices at distance 1
# The number of variables which appear in P and Q optimizations
    def NumVariables(self):
        return len(self.vertices) + len(self.neighbors)

# This is the number of vertices in the support of the prevector
# The number of constraints in the P and Q optimizations
    
    def NumVertexVariables(self):
        return len(self.vertices)

    def setLaplacian(self, v):
	self.Laplacian = v

    def setLaplaceConstraint(self, v):
	self.LaplaceConstraint = v

# Obtain the value of the linearized optimization function

    def ObtainFunctionValueSafe(self,x):
        retval = 0
        for y in x:
            retval += ObjFnSafe(y, self.NumBins)
        return retval
# Obtain the derivative of the linearized optimization function

    def ObtainDerivativeSafe(self, x):
        deriv = np.zeros(self.NumVariables())
        for i in range(self.NumVariables()):
            deriv[i] = DerivSafe(x[i], self.NumBins)
        return deriv

# Obtain the Hessian of the linearized optimization function

    def ObtainHessianSafe(self,x):
        Hessian = np.zeros((self.NumVariables(), self.NumVariables()))
        for i in range(self.NumVariables()):
            Hessian[i,i] = SecondDerivSafe(x[i], self.NumBins)
        return Hessian


# Obtain the objective function (not linearized)        
    def ObtainFunctionValue(self, x):
        retval = self.NumVariables()
        for y in x:
            retval -= np.cos(2 * np.pi * y)
        return retval

# Obtain the derivative (not linearized)
    def ObtainDerivative(self, x):
        deriv = np.zeros(self.NumVariables())
        for i in range(self.NumVariables()):
            deriv[i] = 2 * np.pi * np.sin(2 * np.pi *x[i])
        return deriv

# Obtain the Hessian (not linearized)
    def ObtainHessian(self, x):
        Hessian = np.zeros((self.NumVariables(), self.NumVariables()))
        for i in range(self.NumVariables()):
            Hessian[i,i] = 4 * np.pi**2 * np.cos(2 * np.pi * x[i])
        return Hessian

# These bounds are used for the un-careful P optimization problem which obtain a local, but not necessarily global minimum
    
    def ObtainBounds(self):
        l=[]
        for i in range(self.NumVariables()):
            l.append((0, .5))
        return tuple(l)
       

#These bounds are used for the careful P_j optimization problem, which obtains a global minimum for P_j

    def ObtainBoundsCareful(self, ind, bin):
        lower_bound = .25 + (.25 * bin)/self.NumBins
        upper_bound = .25 + (.25 * (bin+1))/self.NumBins
        l = []
        for i in range(self.NumVariables()):
            if i == ind:
                l.append((lower_bound, upper_bound))
            else:
                l.append((0,.25))
        return tuple(l)

#These bounds are used for the careful Q_j optimization problem, which obtains a global minimum for Q_j    
    def ObtainSignedBoundsCareful(self, ind, bin, eps): 
        lower_bound = .25 + (.25 * bin)/self.NumBins
        upper_bound = .25 + (.25 * (bin+1))/self.NumBins
        l=[]
        if eps < 0:
            temp = - upper_bound
            upper_bound = - lower_bound
            lower_bound = temp
        for i in range(self.NumVariables()):
            if i == ind:
                l.append((lower_bound, upper_bound))
            else:
                l.append((-.25, .25))
        return tuple(l)
    
    def ObtainConstraintMatrix(self):
        nv = self.NumVariables()
        nvv = self.NumVertexVariables()
        returnmat = np.zeros((nvv, nv))
        for i in range(nvv):
            returnmat[i,i] = 6
            for n in self.vertices[i].ObtainNeighbors():
                if not NotPresent(self.vertices, n):
                    ind = bisect.bisect(self.vertices, n)-1
                else:
                    ind = nvv + bisect.bisect(self.neighbors, n)-1
                returnmat[i, ind]=1
        return returnmat
            
# This is the unsafe optimization problem P
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

# This is the safe (binned) optimization problem P_j
    def ObtainValueCareful(self):
        bestValue = self.ObtainOptimizationValue()
        for ind in range(self.NumVariables()):
            for bin in range(self.NumBins):
                c = np.zeros(self.NumVariables())
                c[0] = 1
                b = self.LaplaceConstraint
                A = self.ObtainConstraintMatrix()
                opt_bounds = self.ObtainBoundsCareful(ind,  bin)
                res = linprog(c, A_eq = A, b_eq = b, bounds = opt_bounds)
                if res.status == 0:
                    x0 = res.x
                    def cons(x):
                        return A.dot(x)-b
                    constr = ({'type': 'ineq', 'fun':cons})
                    res = minimize(lambda x: self.ObtainFunctionValueSafe(x), x0, method = 'SLSQP', bounds = opt_bounds, constraints = constr, jac = lambda x: self.ObtainDerivativeSafe(x), hess = lambda x: self.ObtainHessianSafe(x))
                    if res.success:
                        bestValue = min(bestValue, res.fun)
        self.bestValue = bestValue
        return self.bestValue

        

# This gives the bounds for the unsafe signed optimization problem Q

    def ObtainSignedBounds(self):
        l = []
        for i in range(self.NumVariables()):
            l.append((-.5,.5))
        return tuple(l)


    def ObtainSignedConstraintMatrix(self):
        nv = self.NumVariables()
        nvv = self.NumVertexVariables()
        returnmat = np.zeros((nvv, nv))
        for i in range(nvv):
            returnmat[i,i] = 6
            for n in self.vertices[i].ObtainNeighbors():
                if not NotPresent(self.vertices, n):
                    ind = bisect.bisect(self.vertices, n)-1
                else:
                    ind = nvv + bisect.bisect(self.neighbors, n)-1
                returnmat[i, ind]=-1
        return returnmat


# This is the unsafe optimization problem Q

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

# This is the careful signed optimization problem Q_j
        
    def ObtainSignedOptimizationCareful(self):
        bestValue = self.ObtainSignedOptimizationValue()
        for ind in range(self.NumVariables()):
            for eps in [-1,1]:
                for bin in range(self.NumBins):
                    lower_bound = .25 + (.25*bin)/self.NumBins
                    if CosDif(lower_bound) < target_value:
                        opt_bounds = self.ObtainSignedBoundsCareful(ind, bin, eps)
                        b = self.LaplaceConstraint
                        A = self.ObtainSignedConstraintMatrix()
                        c = np.zeros(self.NumVariables())
                        c[0]=1
                        res = linprog(c, A_eq = A, b_eq = b, bounds =opt_bounds)
                        if res.status ==0:
                            x0 = res.x
                            def cons(x):
                                return A.dot(x)-b
                            constr = ({'type': 'eq', 'fun':cons})
                            res = minimize(lambda x: self.ObtainFunctionValueSafe(x), x0, method = 'SLSQP', bounds = opt_bounds, constraints = constr, jac = lambda x: self.ObtainDerivativeSafe(x), hess = lambda x: self.ObtainHessianSafe(x))
                            if res.success:
                                val = res.fun
                                bestValue = min(bestValue, val)
        self.bestValue = bestValue
        return self.bestValue
        
    
            
    def CheckValue(self):
        if self.value < 0:
            return self.ObtainValue()
        else:
            return self.value

# Obtain the ell^2 norm of xi as a Fourier integral

    def ObtainL2(self):
        def FT(x,y):
            rv = 0
            for k in range(self.n):
                rv += self.Laplacian[k] * cmath.exp(-2 *np.pi *1j *(self.vertices[k].x_coor*x + self.vertices[k].y_coor*y))
            return abs(rv)**2/(6 - 2.*(np.cos(2*np.pi *x)+ np.cos(2 * np.pi*y) + np.cos(2*np.pi * (x-y))))**2
        return SingularityIntegral(FT)

# Obtain xi via Fourier integral

    def ObtainXiRepn(self, n):
        def FT(x,y):
            rv = 0
            for k in range(self.n):
                rv += self.Laplacian[k] * cmath.exp(-2 *np.pi *1j *(self.vertices[k].x_coor*x + self.vertices[k].y_coor*y))
            return rv / (6 - 2.*(np.cos(2*np.pi *x)+ np.cos(2 * np.pi*y) + np.cos(2*np.pi * (x-y))))
        time_domain = np.zeros((2*n+1, 2*n+1))
        for a in range(2*n+1):
            for b in range(2*n+1):
                def integrand(x,y):
                    rvalue = FT(x,y) * cmath.exp(-2 * np.pi *1j *((n-a)*x + (n-b)*y))
                    return rvalue.real
                time_domain[a,b] = SingularityIntegral(integrand)
        return time_domain

# THIS FUNCTION RETURNS THE COSINE SUM OVER N_1 V_1 + N_2 V_2 SUCH THAT |N_1|, |N_2| <= N
    
    def ObtainValue(self, n):
        self.value = 0
        xi = self.ObtainXiRepn(n)
        for a in range(2*n+1):
            for b in range(2*n+1):
                self.value += 1 - np.cos(2 * np.pi * xi[a,b])
        return self.value

# THIS FUNCTION RETURNS THE COSINE SUM OVER N_1 V_1 + N_2 V_2 SUCH THAT |N_1|, |N_2| <= N, WITH L_2 CORRECTION.
    
    def ObtainValueL2Correction(self, n):
        self.value = 0
        L2 = 0
        xi = self.ObtainXiRepn(n)
        for a in range(2*n+1):
            for b in range(2*n+1):
                self.value += 1 - np.cos(2 * np.pi * xi[a,b])
                L2 += xi[a,b]**2
        L2_excess = self.ObtainL2() - L2
        self.value += 2*np.pi**2 *L2_excess - np.pi**4 * L2_excess**2/3.
        error = np.pi**4 * L2_excess**2/3.
        return [self.value, error]


# The optimization values in Lemma 13 are obtained here #

#In [16]: c = Configuration([v])

#In [17]: c.setLaplaceConstraint([2])

#In [23]: c.NumBins = 100

#In [24]: c.ObtainValueCareful()
#Out[24]: 1.4322253399783444

#In [25]: c.setLaplaceConstraint([1])

#In [26]: c.ObtainValueCareful()
#Out[26]: 0.4425614825454105

# The optimization values in Lemma 14 are obtained here #

#In [27]: v1 = Vertex(1,1)

#In [28]: v2 = Vertex(2,0)

#In [29]: c1 = Configuration([v,v1])

#In [30]: c2 = Configuration([v,v2])

#In [31]: c1.setLaplaceConstraint([2,1])

#In [32]: c2.setLaplaceConstraint([2,1])

#In [33]: c1.ObtainValueCareful()
#Out[33]: 1.830901686216766

#In [34]: c2.ObtainValueCareful()
#Out[34]: 1.8519792470187433

# The optimization value in Lemma 15 is obtained here #

#In [41]: c = Configuration([v,v1,v2,v3,v4,v5,v6])

#In [42]: c.vertices
#Out[42]: [(-1,0) , (-1,1) , (0,-1) , (0,0) , (0,1) , (1,-1) , (1,0) ]

#In [43]: c.LaplaceConstraint = [1,1,1,2,1,1,1]

#In [44]: c.ObtainValueCareful()
#Out[44]: 1.9233549553724203

# The value in Lemma 16 is obtained here #

#In [45]: c1 = Configuration([v,v1])
#In [48]: c1.LaplaceConstraint = [2,2]

#In [49]: c1.ObtainValueCareful()
#Out[49]: 2.3050944784284324

# The value of [-1,2,-1] is checked here#
#In [15]: c.vertices

#Out[15]: [(-1,0) , (0,0) , (1,0) ]

#In [16]: c.Laplacian = [-1,2,-1]

#In [17]: c.ObtainValue(5)
#Out[17]: 2.231906422129347



# The values in Lemma 17 are obtained here #

# In [50]: v1 = Vertex(1,0)

#In [51]: v2 = Vertex(1,1)

#In [52]: v3 = Vertex(2,0)

#In [53]: c1 = Configuration([v,v1])

#In [54]: c2 = Configuration([v,v2])

#In [55]: c3 = Configuration([v,v3])

#In [56]: c1.LaplaceConstraint = [1,1]

#In [57]: c2.LaplaceConstraint = [1,1]

#In [58]: c3.LaplaceConstraint = [1,1]

#In [59]: c1.ObtainValueCareful()
#Out[59]: 0.6729389504071732

#In [60]: c2.ObtainValueCareful()
#Out[60]: 0.8509744236183668

#In [61]: c3.ObtainValueCareful()
#Out[61]: 0.8677764783693719

# The value in Lemma 18 is obtained here

#In [77]: l
#Out[77]: [(-1,0) , (-1,1) , (0,-1) , (0,0) , (0,1) , (1,-1) , (1,0) ]

#In [78]: c = Configuration(l)

#In [79]: c.LaplaceConstraint = [0,0,0,1,0,0,0]

#In [80]: c.ObtainSignedOptimizationCareful()
#Out[80]: 0.9127672302181694

# The value of the optimization program in Lemma 19 is obtained as follows

#In [81]: v
#Out[81]: (0,0) 

#In [82]: v1 = Vertex(1,0)

#In [83]: c = Configuration([v,v1])

#In [84]: c.LaplaceConstraint = [1,1]

#In [85]: c.ObtainSignedOptimizationCareful()
#Out[85]: 1.1518781502677464

# The value in Lemma 20 is obtained here

#In [89]: l
#Out[89]: 
#[(-1,0) ,
# (-1,1) ,
# (0,-1) ,
# (0,0) ,
# (0,1) ,
# (1,-1) ,
# (1,0) ,
# (1,1) ,
# (2,-1) ,
# (2,0) ]

#In [93]: c.LaplaceConstraint = [0,0,0,1,0,0,-1,0,0,0]

#In [94]: c.ObtainSignedOptimizationCareful()
#Out[94]: 0.9717339887543212


# Two pairs of adjacent vertices are ruled out here:

# In [112]: l = ObtainDist3PairNeighbors()

# In [115]: l1 = Obtain4C2Components(l)

# In [118]: for c in l1:
#     ...:     if c.ObtainValue(5) < 1.7:
#     ...:         print c
#     ...:     

# These routines demonstrate that the best configuration on a single connected component is \nu_0

# In [132]: cc_triangular = RecursivelyBuildComponents() 

# In [147]: len(cc_triangular)
# Out[147]: 8

# In [148]: len(cc_triangular[7])
# Out[148]: 0

# In [133]: cc4 = Obtain4C2Components(cc_triangular[3])

# In [134]: for c in cc4:
#     ...:     val = c.ObtainValue(5)
#     ...:     print val
#     ...:     if val < 1.7:
#     ...:         l.append(c)
#     ...:   

# In [135]: l
# Out[135]: 
# [[[ 1.  1.]
#   [ 1.  1.]]]

# In [139]: l = []

# In [140]: for c in cc6:
#     ...:     val = c.ObtainValue(5)
#     ...:     print val
#     ...:     if val < 1.7:
#     ...:         l.append(c)
#     ...: 

# In [141]: l
# Out[141]: []

# The following rules out configurations with two connected components, one a singleton

# In [151]: cc_triangular_1 = RecursivelyBuildComponents()

# In [152]: len(cc_triangular_1)
# Out[152]: 6

# In [153]: len(cc_triangular_1[4])
# Out[153]: 1

# In [155]: cc_triangular_1[4]
# Out[155]: 
# [[[ 1.  1.  1.]
#   [ 1.  1.  0.]]]

# In [172]: l = ObtainSingletonC2_5(c)

# In [175]: l1 = Obtain6C2Components(l)

# In [177]: for c in l1:
#     ...:     print c.ObtainValue(5)
#     ...:     
#5.47497574651
#8.00858110709
#5.47497574651
#8.00858110709



# The following rules out configurations with two connected components, one a pair of adjacent vertices with opposite sign

# In [161]: cc_triangular_2[3]
# Out[161]: 
# [[[ 1.  1.]
#   [ 1.  1.]]]


# In [164]: len(cc_triangular_2)
# Out[164]: 5

############################

# The value gamma_tri was calculated as follows.

#In [10]: print c.ObtainValueL2Correction(10)
#[1.6941656166936669, 5.294419716854114e-07]







