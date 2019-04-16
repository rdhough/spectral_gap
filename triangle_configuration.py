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

# Cosine difference in optimization functional

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

def SecondDerivSafe(x):
    if (abs(x) < 0.25):
        return 4*np.pi**2 *np.cos(2*np.pi*x)
    else:
        return 0


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

 ## The coordinates give the location of x v_1 + y v_2.  This object can be sorted in a sorted list, and can locate its neighbors in the lattice graph
 
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
    def __str__(self):
        return "(" +str(self.x_coor) + "," +str(self.y_coor) + ") "
    def __repr__(self):
        return "(" +str(self.x_coor) + "," +str(self.y_coor) + ") "


# A configuration consists of a list of vertices in the support of a prevector.  The Laplacian is the function values at the points in the support.  This object can call an optimization program for either Q_j or P_j optimizations, or can use the Fourier integral to obtain the values of $\xi$ in a neighborhood of 0.
    
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
        self.x_max = vertices[0].x_coor
        self.x_min = vertices[0].x_coor
        self.y_max = vertices[0].y_coor
        self.y_min = vertices[0].y_coor
        for ve in vertices:
            self.x_max = max(self.x_max, ve.x_coor)
            self.x_min = min(self.x_min, ve.x_coor)
            self.y_max = max(self.y_max, ve.y_coor)
            self.y_min = min(self.y_min, ve.y_coor)
        self.prevector = np.zeros((self.x_max - self.x_min + 1, self.y_max-self.y_min + 1))
        for ve in vertices:
            self.prevector[ve.x_coor-self.x_min, ve.y_coor-self.y_min] =1
        self.value = -1
        self.bestValue = -1
        self.NumBins = 4

    def __str__(self):
        return str(self.prevector)

    def __repr__(self):
        return str(self.prevector)

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

    def AddElem(self, v):
        newverts = []
        for v_ in self.vertices:
            v1 = Vertex(v_[1], v_[2])
            newverts.append(v1)
        if NotPresent(newverts, v):
	    bisect.insort(newverts, v)
        return Configuration(newverts)


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
    
    def NumVariables(self):
        return len(self.vertices) + len(self.neighbors)

# This is the number of vertices in the support of the prevector
    
    def NumVertexVariables(self):
        return len(self.vertices)

    def setLaplacian(self, v):
	self.Laplacian = v

    def setLaplaceConstraint(self, v):
	self.LaplaceConstraint = v

    def ObtainFunctionValueSafe(self,x):
        retval = 0
        for y in x:
            retval += ObjFnSafe(y, self.NumBins)
        return retval

    def ObtainDerivativeSafe(self, x):
        deriv = np.zeros(self.NumVariables())
        for i in range(self.NumVariables()):
            deriv[i] = DerivSafe(x[i], self.NumBins)
        return deriv

    def ObtainHessianSafe(self,x):
        Hessian = np.zeros((self.NumVariables(), self.NumVariables()))
        for i in range(self.NumVariables()):
            Hessian[i,i] = SecondDerivSafe(x[i], self.NumBins)
        return Hessian


        
    def ObtainFunctionValue(self, x):
        retval = self.NumVariables()
        for y in x:
            retval -= np.cos(2 * np.pi * y)
        return retval

    def ObtainDerivative(self, x):
        deriv = np.zeros(self.NumVariables())
        for i in range(self.NumVariables()):
            deriv[i] = 2 * np.pi * np.sin(2 * np.pi *x[i])
        return deriv

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

    def ObtainL2(self):
        def FT(x,y):
            rv = 0
            for k in range(self.n):
                rv += self.Laplacian[k] * cmath.exp(-2 *np.pi *1j *(self.vertices[k].x_coor*x + self.vertices[k].y_coor*y))
            return abs(rv)**2/(6 - 2.*(np.cos(2*np.pi *x)+ np.cos(2 * np.pi*y) + np.cos(2*np.pi * (x-y))))**2
        return SingularityIntegral(FT)
    #dblquad(lambda t,x: FT(t,x), 0, 1, lambda x:0, lambda x:1, epsabs = 1.0e-11)[0]
    
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
                #integrate.nquad( lambda x, t: integrand(x,t), [[0,1],[0,1]], full_output = True)
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



