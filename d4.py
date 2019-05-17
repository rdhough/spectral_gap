##################################################################################################
# This file provides the D4 lattice configuration class, which is used to calculate the extremal #
# configurations which give the spectral factors of the D4 lattice. These functions are used in  #
# the iPython environment.  Example usage which obtains the spectral factors is given at the end #
# of the file.                                                                                   #
##################################################################################################

# Import standard libraries


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

########################################################################################################
# The next few functions provide a method of integrating in polar coordinates to avoid the singularity #
# of the Green's function at 0. Since the Fourier integrals are over (R/Z)^4, a spherical neighborhood #
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

def SphericalCoordinateIntegrate(integrand):
    def spherical_integrand(r,thet1, thet2, thet3):
        if r == 0:
            return 0
        else:
            return integrand(r*np.cos(thet1), r*np.sin(thet1)*np.cos(thet2),r*np.sin(thet1)*np.sin(thet2)*np.cos(thet3), r*np.sin(thet1)*np.sin(thet2)*np.sin(thet3) )*r**3 * np.sin(thet1)**2 *np.sin(thet2)
        
    return integrate.nquad(lambda r, thet1, thet2, thet3: spherical_integrand(r,thet1, thet2, thet3), [[0, .5],[0,np.pi],[0, np.pi], [0, 2*np.pi]])[0]

# Integration over the torus (R/Z)^4 with singularity at 0

def SingularityIntegral(integrand):
    def smooth_integrand(x,y,z,w):
        if x**2 + y**2+z**2+w**2 == 0:
            return 0
        else:
            return integrand(x,y,z,w)* (1-BumpFunction(4*(x**2 + y**2+z**2+w**2)))
    def singular_integrand(x,y,z,w):
        return integrand(x,y,z,w) * BumpFunction(4*(x**2 + y**2+ z**2+w**2))
    SmoothIntegral = integrate.nquad(lambda x,y,z,w: smooth_integrand(x,y,z,w), [[-.5,.5],[-.5,.5], [-.5,.5],[-.5,.5]])[0]
    SingularIntegral = SphericalCoordinateIntegrate(singular_integrand)
    return SmoothIntegral + SingularIntegral

# The following function is used in defining f(xi)

def CosDif(x):
    return 1 - np.cos(2 *np.pi * x)

#######################################
# Helper functions to work with lists #
#######################################

def SortList(l):
    returnlist = []
    for e in l:
        bisect.insort(returnlist,e)
    return returnlist

def NotPresent(sort_list, elem):
    index = bisect.bisect(sort_list, elem)
    if index == 0 or elem != sort_list[index-1]:
        return True
    return False

def MergeList(l1, l2):
    returnlist = SortList(l1)
    for e in l2:
        if NotPresent(returnlist, e):
	    bisect.insort(returnlist,e)
    return returnlist


# This implements a point in the D4 lattice

class Vertex(object):
    def __init__(self, w_coor, x_coor, y_coor, z_coor):
	self.w_coor = w_coor
        self.x_coor = x_coor
        self.y_coor = y_coor
        self.z_coor = z_coor
    def __lt__(self, other):
	if self.w_coor < other.w_coor:
	    return True
        elif self.w_coor == other.w_coor and self.x_coor < other.x_coor:
            return True
        elif self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor < other.y_coor:
            return True
        elif self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor < other.z_coor:
            return True
        else:
            return False
    def __le__(self, other):
	if self.w_coor < other.w_coor:
	    return True
        elif self.w_coor == other.w_coor and self.x_coor < other.x_coor:
            return True
        elif self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor < other.y_coor:
            return True
        elif self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor <= other.z_coor:
            return True
        else:
            return False
    def __gt__(self, other):
	if self.w_coor > other.w_coor:
	    return True
        if self.w_coor == other.w_coor and self.x_coor > other.x_coor:
            return True
        elif self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor > other.y_coor:
            return True
        elif self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor > other.z_coor:
            return True
        else:
            return False
    def __ge__(self, other):
	if self.w_coor > other.w_coor:
	    return True
        elif self.w_coor == other.w_coor and self.x_coor > other.x_coor:
            return True
        elif self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor > other.y_coor:
            return True
        elif self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor >= other.z_coor:
            return True
        else:
            return False
    def __eq__(self, other):
        if self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor == other.z_coor:
            return True
        else:
            return False
    def __ne__(self, other):
        if self.w_coor == other.w_coor and self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor == other.z_coor:
            return False
        else:
            return True
    def __getitem__(self, key):
        if key == 1:
            return self.w_coor
        if key == 2:
            return self.x_coor
        if key == 3:
            return self.y_coor
	if key == 4:
	    return self.z_coor
# Obtains the 24 nearest neighbors of the point        
    def ObtainNeighbors(self):
	v1 = Vertex(self.w_coor + 1, self.x_coor, self.y_coor, self.z_coor)
	v2 = Vertex(self.w_coor - 1, self.x_coor, self.y_coor, self.z_coor)
	v3 = Vertex(self.w_coor, self.x_coor+1, self.y_coor, self.z_coor)
	v4 = Vertex(self.w_coor, self.x_coor-1, self.y_coor, self.z_coor)
	v5 = Vertex(self.w_coor, self.x_coor, self.y_coor +1, self.z_coor)
	v6 = Vertex(self.w_coor, self.x_coor, self.y_coor-1, self.z_coor)
	v7 = Vertex(self.w_coor, self.x_coor, self.y_coor, self.z_coor+1)
	v8 = Vertex(self.w_coor, self.x_coor, self.y_coor, self.z_coor-1)
	v9 = Vertex(self.w_coor+.5, self.x_coor +.5, self.y_coor + .5, self.z_coor + .5)
	v10 = Vertex(self.w_coor+.5, self.x_coor +.5, self.y_coor + .5, self.z_coor - .5)
	v11 = Vertex(self.w_coor+.5, self.x_coor +.5, self.y_coor - .5, self.z_coor + .5)
	v12 = Vertex(self.w_coor+.5, self.x_coor +.5, self.y_coor - .5, self.z_coor - .5)
	v13 = Vertex(self.w_coor+.5, self.x_coor -.5, self.y_coor + .5, self.z_coor + .5)
	v14 = Vertex(self.w_coor+.5, self.x_coor -.5, self.y_coor + .5, self.z_coor - .5)
	v15 = Vertex(self.w_coor+.5, self.x_coor -.5, self.y_coor - .5, self.z_coor + .5)
	v16 = Vertex(self.w_coor+.5, self.x_coor -.5, self.y_coor - .5, self.z_coor - .5)
	v17 = Vertex(self.w_coor-.5, self.x_coor +.5, self.y_coor + .5, self.z_coor + .5)
	v18 = Vertex(self.w_coor-.5, self.x_coor +.5, self.y_coor + .5, self.z_coor - .5)
	v19 = Vertex(self.w_coor-.5, self.x_coor +.5, self.y_coor - .5, self.z_coor + .5)
	v20 = Vertex(self.w_coor-.5, self.x_coor +.5, self.y_coor - .5, self.z_coor - .5)
	v21 = Vertex(self.w_coor-.5, self.x_coor -.5, self.y_coor + .5, self.z_coor + .5)
	v22 = Vertex(self.w_coor-.5, self.x_coor -.5, self.y_coor + .5, self.z_coor - .5)
	v23 = Vertex(self.w_coor-.5, self.x_coor -.5, self.y_coor - .5, self.z_coor + .5)
	v24 = Vertex(self.w_coor-.5, self.x_coor -.5, self.y_coor - .5, self.z_coor - .5)	
        return SortList([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20,v21,v22,v23,v24])

# Obtains the graph-distance 2 neighborhood of the point    
    def ObtainDist2Neighbors(self):
	l = self.ObtainNeighbors()
        returnlist = l
	for e in l:
	    returnlist = MergeList(returnlist, e.ObtainNeighbors())
        return returnlist

# Obtains the graph-distance 3 neighborhood of the point    
    def ObtainDist3Neighbors(self):
        l = self.ObtainDist2Neighbors()
        returnlist = l
        for e in l:
            returnlist = MergeList(returnlist, e.ObtainNeighbors())
        return returnlist

# Obtains the graph-distance 4 neighborhood of the point    
    def ObtainDist4Neighbors(self):
        l = self.ObtainDist3Neighbors()
        returnlist = l
        for e in l:
            returnlist = MergeList(returnlist, e.ObtainNeighbors())
        return returnlist
    
# Display
    def __str__(self):
        return "(" +str(self.w_coor)+","+str(self.x_coor) + "," +str(self.y_coor) +","+ str(self.z_coor)+") "
    def __repr__(self):
        return "(" +str(self.w_coor)+","+str(self.x_coor) + "," +str(self.y_coor) +","+ str(self.z_coor)+") "

    
# Multiplication on the right by quaternion q
    def right_mult(self, q):
        new_w = self.w_coor * q.w_coor - self.x_coor * q.x_coor - self.y_coor * q.y_coor - self.z_coor *q.z_coor
        new_x = self.w_coor * q.x_coor +self.x_coor *q.w_coor + self.y_coor*q.z_coor - self.z_coor * q.y_coor
        new_y = self.w_coor * q.y_coor + self.y_coor*q.w_coor + self.z_coor *q.x_coor - self.x_coor *q.z_coor
        new_z = self.w_coor * q.z_coor + self.z_coor *q.w_coor + self.x_coor *q.y_coor - self.y_coor *q.x_coor
        return Vertex(new_w, new_x, new_y, new_z)

# Multiplies the vector by each quaternion unit, then returns the least such element    
    def Normalize(self):
        v0 = Vertex(0,0,0,0)
        l = v0.ObtainNeighbors()
        lis = []
        for v in l:
            lis.append(self.right_mult(v))
        l2 = SortList(lis)
        v1 = l2[0]
        v2 = Vertex(abs(v1.w_coor), abs(v1.x_coor), abs(v1.y_coor), abs(v1.z_coor))
        return v2
    

#######################################################################################################
# This class implements configurations of vertices in the D4 lattice.  A configuration consists of a  #
# list of points in the lattice, together with their "Laplacian" prevector, which is an integer       #
# valued function on the vertices.                                                                    #
#######################################################################################################
   
    
class Configuration(object):
    def __init__(self, vertices):
        self.vertices = SortList(vertices)
        self.n = len(vertices)
# The Laplacian (prevector values) and Laplace constraint are intialized to 0 and must be set prior to use
	self.LaplaceConstraint = np.ones(self.n)
	self.Laplacian = np.zeros(self.n)
        self.neighbors = []
        self.dist2neighbors = []
# Set the neighbors of the point
        for v in vertices:
            for n in v.ObtainNeighbors():
                if NotPresent(self.neighbors, n) and NotPresent(self.vertices, n):
                    bisect.insort(self.neighbors, n)
        for v in vertices:
            for n in v.ObtainDist2Neighbors():
                if NotPresent(self.dist2neighbors, n) and NotPresent(self.vertices, n):
                    bisect.insort(self.dist2neighbors, n)
# Bounds for the values of the coordinates of the configuration
	self.w_max = vertices[0].w_coor
	self.w_min = vertices[0].w_coor
        self.x_max = vertices[0].x_coor
        self.x_min = vertices[0].x_coor
        self.y_max = vertices[0].y_coor
        self.y_min = vertices[0].y_coor
        self.z_max = vertices[0].z_coor
        self.z_min = vertices[0].z_coor
        for ve in vertices:
	    self.w_max = max(self.w_max, ve.w_coor)
	    self.w_min = min(self.w_min, ve.w_coor)
            self.x_max = max(self.x_max, ve.x_coor)
            self.x_min = min(self.x_min, ve.x_coor)
            self.y_max = max(self.y_max, ve.y_coor)
            self.y_min = min(self.y_min, ve.y_coor)
            self.z_max = max(self.z_max, ve.z_coor)
            self.z_min = min(self.z_min, ve.z_coor)
        self.value = -1
        self.NumBins = 5
        v0 = Vertex(0,0,0,0)
        v0_neighbors = v0.ObtainNeighbors()
   
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
        if (self.ge(other) and self.le(other)):
            return True
        else:
            return False

    def __ne__(self, other):
        if (self.ge(other) and self.le(other)):
            return False
        else:
            return True

# The number of variables which appear in the programs P and Q
    def NumVariables(self):
        return len(self.vertices) + len(self.neighbors)

# The number of constraints which appear in the programs P and Q
    def NumVertexVariables(self):
        return len(self.vertices)

# Set the values of the prevector
    def setLaplacian(self, v):
	self.Laplacian = v

# Set the constraint for the optimization programs
    def setLaplaceConstraint(self, v):
	self.LaplaceConstraint = v

# Sum of squares#
    def ObtainFunctionValue(self, x):
        retval = 0
        for y in x:
            retval += y*y
        return retval

# Gradient for sum of squares#
    def ObtainDerivative(self, x):
        deriv = np.zeros(self.NumVariables())
        for i in range(self.NumVariables()):
            deriv[i] = 2 *x[i]
        return deriv

# Hessian for sum of squares#
    def ObtainHessian(self, x):
        Hessian = np.zeros((self.NumVariables(), self.NumVariables()))
        for i in range(self.NumVariables()):
            Hessian[i,i] = 2 
        return Hessian

# Bounds for the program Q'#
    def ObtainBounds(self):
        l=[]
        for i in range(self.NumVariables()):
            l.append((-.5, .5))
        return tuple(l)
   
# Constraint Matrix for the program Q' #
    def ObtainConstraintMatrix(self):
        nv = self.NumVariables()
        nvv = self.NumVertexVariables()
        returnmat = np.zeros((nvv, nv))
        for i in range(nvv):
            returnmat[i,i] = 24
            for n in self.vertices[i].ObtainNeighbors():
                if not NotPresent(self.vertices, n):
                    ind = bisect.bisect(self.vertices, n)-1
                else:
                    ind = nvv + bisect.bisect(self.neighbors, n)-1
                returnmat[i, ind]=-1
        return returnmat
            
# This optimization program is the program Q'#
    def ObtainOptimizationValue(self):
        x0 = np.zeros(self.NumVariables())
        for j in range(self.NumVariables()):
            x0[j] = .125
        A = self.ObtainConstraintMatrix()
        def cons(x):
            return A.dot(x) - self.LaplaceConstraint
        constr = ({'type': 'eq', 'fun': cons})
        res = minimize(lambda x: self.ObtainFunctionValue(x),x0, method ='SLSQP', bounds = self.ObtainBounds(), constraints =constr, jac = lambda x: self.ObtainDerivative(x), hess = lambda x: self.ObtainHessian(x) )
        if res.success:
            self.value = res.fun
            return   res.fun
        else:
            return -1


#########################################################################################################################
# This calculation checks that for each i= 1, 2, 3, 4 there is a single element of norm i up to multiplication by a unit.  The first and second Lemma numbers corresponding to the corresponding lemma in the short or long version of the paper.  The second number is the lemma number in the version in the repository.#
#########################################################################################################################

################################################
# The values in Lemma 15/33 are confirmed here #
################################################

#In [50]: l = v0.ObtainNeighbors()

#In [51]: bisect.insort(l,v0)

#In [52]: len(l)
#Out[52]: 25

#In [53]: c = Configuration(l)

#In [54]: c.LaplaceConstraint = np.zeros(25)

#In [55]: c.LaplaceConstraint[12] = 1

#In [56]: c.ObtainOptimizationValue()
#Out[56]: 0.002063197026022368


#In [57]: l = v0.ObtainDist2Neighbors()

#In [61]: len(l)
#Out[61]: 169

#In [62]: l[84]
#Out[62]: (0,0,0,0) 

#In [63]: c = Configuration(l)

#In [64]: c.LaplaceConstraint = np.zeros(169)

#In [65]: c.LaplaceConstraint[84] = 1

#In [66]: c.ObtainOptimizationValue()
#Out[66]: 0.0023398440686087842

############################################
#The table in Lemma 15/33 is confirmed here#
############################################
#In [67]: v0
#Out[67]: (0,0,0,0) 

#In [68]: v1 = Vertex(1,0,0,0)

#In [69]: v2 = Vertex(1,1,0,0)

#In [70]: v3 = Vertex(1,1,1,0)

#In [71]: v4 = Vertex(2,0,0,0)

#In [72]: c1 = Configuration([v0,v1])

#In [73]: c2 = Configuration([v0,v2])

#In [74]: c3 = Configuration([v0,v3])

#In [75]: c4 = Configuration([v0,v4])

### The first column in the table is obtained here ###

#In [76]: c1.ObtainOptimizationValue()
#Out[76]: 0.0035714285714691296

#In [77]: c2.ObtainOptimizationValue()
#Out[77]: 0.003300330033025527

#In [78]: c3.ObtainOptimizationValue()
#Out[78]: 0.003322259136230209

#In [79]: c4.ObtainOptimizationValue()
#Out[79]: 0.0033277870216472487

#In [80]: c1.LaplaceConstraint = [1,-1]

#In [81]: c2.LaplaceConstraint = [1,-1]

#In [82]: c3.LaplaceConstraint = [1,-1]

#In [83]: c4.LaplaceConstraint = [1,-1]

### The second column in the table is obtained here ###

#In [84]: c1.ObtainOptimizationValue()
#Out[84]: 0.003125000000023783

#In [85]: c2.ObtainOptimizationValue()
#Out[85]: 0.00336700336702744

#In [86]: c3.ObtainOptimizationValue()
#Out[86]: 0.0033444816053692303

#In [87]: c4.ObtainOptimizationValue()
#Out[87]: 0.0033388981636228786


##########################################
# The case of gamma_D4,0 is checked here #
##########################################

#def D4_integrand(w,x,y,z):
#    return 2*(1-np.cos(2*np.pi * w))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4 = integrate.nquad(lambda w,x,y,z: D4_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [39]: L2_D4
#Out[39]: (0.003839735406071655, 1.7839615201749152e-08)

#In [170]: xi2 = 0.00383973

# Calculation of gamma_D4,0

#In [73]: 2*np.pi**2 * alpha - np.pi**4*alpha**2/6.
#Out[73]: 0.07555387329526285

#In [74]: np.pi**4*alpha**2/6.
#Out[74]: 0.00023935891872682344

# Calculation of Gamma_D4,0

#In [75]: 4/(.075554-.00024)
#Out[75]: 53.110975383062915

#In [76]: 4/(.075554+.00024)
#Out[76]: 52.77462595983851

#In [77]: (53.110975383062915+52.77462595983851 )/2.
#Out[77]: 52.94280067145071

#In [78]: (53.110975383062915-52.77462595983851 )/2.
#Out[78]: 0.16817471161220254

####################################################################################
# The following calculation shows that no other pair of nodes gives a better value #
####################################################################################

#
#In [9]: v0 = Vertex(0,0,0,0)

#In [10]: l = v0.ObtainDist4Neighbors()
#In [11]: l1 = []


#In [12]: for v in l:
#    ...:     v1 = v.Normalize()
#    ...:     if NotPresent(l1, v1):
#    ...:         bisect.insort(l1,v1)
#    ...:         

### l1 is a list of vertices whose distance from 0 is at most 4, taken up to multiplication by a unit

#In [13]: len(l1)
#Out[13]: 25

#In [14]: l1
#Out[14]: 
#[(0,0,0,0) ,
# (1,0,0,0) ,
# (1,1,0,0) ,
# (1.5,0.5,0.5,0.5) ,
# (2,0,0,0) ,
# (2,0,0,1) ,
# (2,0,1,0) ,
# (2,1,0,0) ,
# (2.0,1.0,0.0,1.0) ,
# (2.0,1.0,1.0,0.0) ,
# (2,2,0,0) ,
# (2.5,0.5,0.5,0.5) ,
# (2.5,0.5,0.5,1.5) ,
# (2.5,0.5,1.5,0.5) ,
# (2.5,1.5,0.5,0.5) ,
# (3,0,0,0) ,
# (3,0,0,1) ,
# (3,0,1,0) ,
# (3.0,0.0,1.0,1.0) ,
# (3,1,0,0) ,
# (3.0,1.0,0.0,1.0) ,
# (3.0,1.0,1.0,0.0) ,
# (3.0,1.0,1.0,1.0) ,
# (3.5,0.5,0.5,0.5) ,
# (4,0,0,0) ]


#In [17]: for i in range(24):
#    ...:     ...:     v1 = l1[i+1]
#    ...:     ...:     l4 = v1.ObtainDist2Neighbors()
#    ...:     ...:     l5 = MergeList(l2, l4)
#    ...:     ...:     i1 = bisect.bisect(l5,v0)-1
#    ...:     ...:     i2 = bisect.bisect(l5, v1)-1
#    ...:     ...:     c = Configuration(l5)
#    ...:     ...:     c.setLaplaceConstraint(np.zeros(len(l5)))
#    ...:     ...:     c.LaplaceConstraint[i1]=1
#    ...:     ...:     c.LaplaceConstraint[i2]=-1
#    ...:     ...:     val = c.ObtainOptimizationValue()
#    ...:     ...:     print i
#    ...:     ...:     print val

#0
#0.00366984749863
#1
#0.0041780902205
#2
#0.00436565218923
#3
#0.00447307451018
#4
#0.00455399436812
#5
#0.00455399436812
#6
#0.00455399436812
#7
#0.00459878797526
#8
#0.00459878797526
#9
#0.00465353773964
#10
#0.00462935396815
#11
#0.00466502073218
#12
#0.00466502073218
#13
#0.00466502073218
#14
#0.00466771619724
#15
#0.00467523438077
#16
#0.00467523438077
#17
#0.0046816882096
#18
#0.00467523438077
#19
#0.0046816882096
#20
#0.0046816882096
#21
#0.00468559301053
#22
#0.00468707857354
#23
#0.00468603244899

### The best value corresponds to (1,0,0,0), and all other values exceed alpha, so are not extremal. ###

############################################
# The case of gamma_D4,1 is checked here   #
############################################


#def D4_1_reflection_integrand(w,x,y,z):
#    return (1-np.cos(2*np.pi * (w+x)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_1_reflections = integrate.nquad(lambda w,x,y,z: D4_1_reflection_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [41]: L2_D4_1_reflections
#Out[41]: (0.002242187555185342, 1.5801954702554342e-08)


#In [173]: xi2 = 0.0022421875

#Calculation of gamma_D4,1

#In [174]: 2*np.pi**2 * xi2 - np.pi**4*xi2**2/3.
#Out[174]: 0.04409576892600773

#In [175]: np.pi**4 *xi2**2/3.
#Out[175]: 0.0001632383101273548

#In [182]: 3/(.0440957-.000163239)
#Out[182]: 68.28663661705635

#In [183]: 3/(.0440957+.000163239)
#Out[183]: 67.78291725429749

#Calculation of Gamma_D4,1

#In [184]: (68.28663661705635 +67.78291725429749)/2
#Out[184]: 68.03477693567692

#In [185]: (68.28663661705635 -67.78291725429749)/2
#Out[185]: 0.25185968137942893

###########################################################################################################################
# The claimed value is the best one by using the program from gamma_D4,0 with (0,0,0,0) and (2,2,0,0), which has value 0.00465353773964 which is more than 2 x alpha.#
###########################################

##########################################
# The case of gamma_D4,2 is checked here #
##########################################

#def D4_2_reflections_integrand(w,x,y,z):
#    return (1-np.cos(2*np.pi * (w+x)))*(1-np.cos(2*np.pi*(w-x)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_2_reflections = integrate.nquad(lambda w,x,y,z: D4_2_reflections_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [12]: L2_D4_2_reflections
#Out[12]: (0.0019800331698782685, 1.4899895256774435e-08)


#In [186]: xi2 = 0.00198003

# Calculation of gamma_D4,2

#In [187]: 2*np.pi**2 * xi2 - np.pi**4 * xi2**2/3.
#Out[187]: 0.03895692754698546

#In [188]: np.pi**4 * xi2**2/3.
#Out[188]: 0.000127298057592462

#In [189]: 2./(.0389569275469 - .0001272980575)
#Out[189]: 51.507058560679155

#In [190]: 2./(.0389569275469 + .0001272980575)
#Out[190]: 51.17153964475236

#Calculation of Gamma_D4,2

#In [191]: (51.5070585606 + 51.1715396447)/2
#Out[191]: 51.33929910265

#In [192]: (51.5070585606 - 51.1715396447)/2
#Out[192]: 0.1677594579500017


################################################################################################################################
#The alternative configuration has distance 1 from one hyperplane and distance 2 from the other.  This gives the value below.###
################################################################################################################################

#def D4_2_reflections_double_integrand(w,x,y,z):
#    return (1-np.cos(2*np.pi * (w+x)))*(1-np.cos(4*np.pi*(w-x)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_2_reflections_double = integrate.nquad(lambda w,x,y,z: D4_2_reflections_double_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [128]: L2_D4_2_reflections_double
#Out[128]: (0.0021612824492875534, 1.4899617198582899e-08)




##########################################
# The case of gamma_D4,3 is checked here #
##########################################


#def D4_3_reflections_integrand(w,x,y,z):
#    return (1-np.cos(2*np.pi * (w+x)))*(1-np.cos(2*np.pi*(w-x)))*(1-np.cos(2*np.pi*(2*z-w-x)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_3_reflections = integrate.nquad(lambda w,x,y,z: D4_3_reflections_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [5]: L2_D4_3_reflections
#Out[5]: (0.001873799409647135, 1.4899291519171603e-08)

#################################################
############### This is the best value ##########
#################################################

##################################
#The calculation of alpha is here#
##################################

#In[1]:= Solve[x*(1 - Pi^2/3 *x) - 0.0018737 == 0, x] 

#Out[1]= {{x -> 0.00188539}, {x -> 0.302078}}


#In [193]: xi2 = 0.0018738

# Calculation of gamma_D4,3

#In [194]: 2*np.pi**2 * xi2 - np.pi**4 * xi2**2/3.
#Out[194]: 0.036873324241847194

#In [195]: np.pi**4 * xi2**2/3.
#Out[195]: 0.00011400521167528427

#In [196]: 1./(.0368733242418-.0001140052116)
#Out[196]: 27.20398599273397

#In [197]: 1./(.0368733242418+.0001140052116)
#Out[197]: 27.036285527450442

# Calculation of Gamma_D4,3

#In [198]: (27.20398599273397 + 27.036285527450442)/2
#Out[198]: 27.120135760092204

#In [199]: (27.20398599273397 - 27.036285527450442)/2
#Out[199]: 0.08385023264176361

################################################################################################
# The following two programs show that the optimum has distance at most 2 from all hyperplanes #
################################################################################################


#In [377]: v0 = Vertex(0,0,0,0)

#In [378]: v1 = Vertex(1,1,0,0)

#In [379]: v2 = Vertex(1,-1,0,0)

#In [380]: v3 = Vertex(2,0,0,0)

#In [381]: l = v0.ObtainDist2Neighbors()

#In [382]: l = MergeList(l, v1.ObtainDist2Neighbors())

#In [383]: l = MergeList(l, v2.ObtainDist2Neighbors())

#In [384]: l = MergeList(l, v3.ObtainDist2Neighbors())

#In [385]: i0 = bisect.bisect_left(l, v0)

#In [386]: i1 = bisect.bisect_left(l, v1)

#In [387]: i2 = bisect.bisect_left(l, v2)

#In [388]: i3 = bisect.bisect_left(l, v3)

#In [389]: l[i0]
#Out[389]: (0,0,0,0) 

#In [390]: l[i1]
#Out[390]: (1,1,0,0) 

#In [391]: l[i2]
#Out[391]: (1,-1,0,0) 

#In [392]: l[i3]
#Out[392]: (2,0,0,0) 

#In [393]: c = Configuration(l)

#In [394]: c.LaplaceConstraint = np.zeros(len(l))

#In [395]: c.LaplaceConstraint[i0]=1

#In [396]: c.LaplaceConstraint[i1]=-1

#In [397]: c.LaplaceConstraint[i2]=-1

#In [398]: c.LaplaceConstraint[i3]=1

#In [399]: c.ObtainOptimizationValue()
#Out[399]: 0.007793588373753282

#In [400]: 0.007793588373753282/4
#     ...: 
#     ...: 
#Out[400]: 0.0019483970934383206

# This value exceeds alpha, which rules out (0,0,0,0), (1,1,0,0), (1,-1,0,0), (2,0,0,0)

##############################################
# The second optimization program starts here#
##############################################

#In [401]: v0
#Out[401]: (0,0,0,0) 

#In [402]: v1
#Out[402]: (1,1,0,0) 

#In [403]: v2 = Vertex(2,-2,0,0)

#In [404]: v3 = Vertex(3,-1,0,0)

#In [405]: l = v0.ObtainDist2Neighbors()

#In [406]: l = MergeList(l, v1.ObtainDist2Neighbors())

#In [407]: l = MergeList(l, v2.ObtainDist2Neighbors())

#In [408]: l = MergeList(l, v3.ObtainDist2Neighbors())

#In [409]: i0 = bisect.bisect_left(l, v0)

#In [410]: i1 = bisect.bisect_left(l, v1)

#In [411]: i2 = bisect.bisect_left(l, v2)

#In [412]: i3 = bisect.bisect_left(l, v3)

#In [413]: c = Configuration(l)

#In [414]: c.LaplaceConstraint = np.zeros(len(l))

#In [415]: c.LaplaceConstraint[i0] = 1

#In [416]: c.LaplaceConstraint[i1] = -1

#In [418]: c.LaplaceConstraint[i2] = -1

#In [419]: c.LaplaceConstraint[i3] = 1

#In [420]: c.ObtainOptimizationValue()
#Out[420]: 0.008328160850749503

#In [421]: 0.008328160850749503/4
#Out[421]: 0.002082040212687376

# This value exceeds alpha, which rules out (0,0,0,0), (1,1,0,0), (2,-2,0,0), (3,-1,0,0)

#####################################################################################################################################
#The remaining possibilities have one or two distance 2 displacements from the reflecting hyperplanes.  These cases are checked here#
#####################################################################################################################################

#def D4_3_reflections_double_integrand(w,x,y,z):
#    return (1-np.cos(2*np.pi * (w+x)))*(1-np.cos(4*np.pi*(w-x)))*(1-np.cos(2*np.pi*(2*z-w-x)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_3_reflections_double = integrate.nquad(lambda w,x,y,z: D4_3_reflections_double_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [130]: L2_D4_3_reflections_double
#Out[130]: (0.0019652726551285865, 1.4896357740968197e-08)

#def D4_3_reflections_double_double_integrand(w,x,y,z):
#    return (1-np.cos(2*np.pi * (w+x)))*(1-np.cos(4*np.pi*(w-x)))*(1-np.cos(4*np.pi*(2*z-w-x)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_3_reflections_double_double = integrate.nquad(lambda w,x,y,z: D4_3_reflections_double_double_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [110]: L2_D4_3_reflections_double_double
#Out[110]: (0.0021227270663249516, 1.4898837717854888e-08)



##############################################################################################
# This calculation determines the best configuration for 4 reflecting hyperplanes, gamma_D4,4#
##############################################################################################

#def D4_4_reflections_integrand(w,x,y,z):
#    return (1-np.cos(2*np.pi * (w+x)))*(1-np.cos(2*np.pi*(w-x)))*(1-np.cos(2*np.pi*(2*z-w-x)))*(1-np.cos(2*np.pi*(2*z-w-x-2*y)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_4_reflections = integrate.nquad(lambda w,x,y,z: D4_4_reflections_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [3]: L2_D4_4_reflections
#Out[3]: (0.0018170734649653364, 1.4879547893881528e-08)

##############################
# This value is the best one #
##############################

# Determination of alpha

#In[2]:= Solve[x*(1 - Pi^2/3 *x) - 0.0018171 == 0, x] 


#Out[2]= {{x -> 0.00182809}, {x -> 0.302135}}

# alpha = 0.00182809

#In [200]: xi2 = 0.001817073

# Calculation of gamma_D4,4

#In [201]: 2*np.pi**2 * xi2 - np.pi**4 * xi2**2/3.
#Out[201]: 0.035760376394485836

#In [202]: np.pi**4 * xi2**2/3.
#Out[202]: 0.00010720696131544611

################################################################################
# Next rule out other configurations at distance at most 2 from each hyperplane#
################################################################################
####################
# One distance of 2#
####################

#def D4_4_reflections_double_integrand(w,x,y,z):
#    return (1-np.cos(4*np.pi * (w+x)))*(1-np.cos(2*np.pi*(w-x)))*(1-np.cos(2*np.pi*(2*z-w-x)))*(1-np.cos(2*np.pi*(2*z-w-x-2*y)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_4_reflections_double = integrate.nquad(lambda w,x,y,z: D4_4_reflections_double_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [7]: L2_D4_4_reflections_double
#Out[7]: (0.0018694703725885338, 1.4890326042522233e-08)

# This exceeds alpha = 0.00182809

#####################
# Two distances of 2#
#####################

#def D4_4_reflections_double_double_integrand(w,x,y,z):
#    return (1-np.cos(4*np.pi * (w+x)))*(1-np.cos(4*np.pi*(w-x)))*(1-np.cos(2*np.pi*(2*z-w-x)))*(1-np.cos(2*np.pi*(2*z-w-x-2*y)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_4_reflections_double_double = integrate.nquad(lambda w,x,y,z: D4_4_reflections_double_double_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [5]: L2_D4_4_reflections_double_double
#Out[5]: (0.001955095910976028, 1.4899986471962798e-08)

# This exceeds alpha = 0.00182809

########################
# Three distances of 2 #
########################

#def D4_4_reflections_double_double_double_integrand(w,x,y,z):
#    return (1-np.cos(4*np.pi * (w+x)))*(1-np.cos(4*np.pi*(w-x)))*(1-np.cos(4*np.pi*(2*z-w-x)))*(1-np.cos(2*np.pi*(2*z-w-x-2*y)))/(24 - 2*(np.cos(2*np.pi * w) + np.cos(2*np.pi *x) + np.cos(2*np.pi * y) + np.cos(2*np.pi * (2*z - w-x-y)) + np.cos(2*np.pi * (z)) + np.cos(2*np.pi * (z-w-x-y)) + np.cos(2*np.pi *(z-w)) + np.cos(2*np.pi * (z-x)) + np.cos(2*np.pi * (z-y)) + np.cos(2*np.pi * (z-w-x)) + np.cos(2*np.pi * (z-w-y)) + np.cos(2*np.pi * (z-x-y))))**2

#L2_D4_4_reflections_double_double_double = integrate.nquad(lambda w,x,y,z: D4_4_reflections_double_double_double_integrand(w,x,y,z), [[0,1],[0,1],[0,1],[0,1]])

#In [112]: L2_D4_4_reflections_double_double_double
#Out[112]: (0.002097967486433312, 1.4899806817847837e-08)

# This exceeds alpha = 0.00182809

#######################################################################################
##### These calculations show that distance 2 to each hyperplane suffices #############
#######################################################################################

#####################################################################
# The first case considered is a node with distance 1 to P1, P2, P3 #
#####################################################################


#In [422]: v0 = Vertex(0,0,0,0)

#In [423]: v1 = Vertex(1,1,0,0)

#In [424]: v2 = Vertex(1,-1,0,0)

#In [425]: v3 = Vertex(2,0,0,0)

#In [426]: v4 = Vertex(0,0,1,1)

#In [427]: v5 = Vertex(1,1,1,1)

#In [428]: v6 = Vertex(1,-1,1,1)

#In [429]: v7 = Vertex(2,0,1,1)

#In [430]: l = v0.ObtainDist2Neighbors()

#In [431]: l = MergeList(l, v1.ObtainDist2Neighbors())

#In [432]: l = MergeList(l, v2.ObtainDist2Neighbors())

#In [433]: l = MergeList(l, v3.ObtainDist2Neighbors())

#In [434]: l = MergeList(l, v4.ObtainDist2Neighbors())

#In [435]: l = MergeList(l, v5.ObtainDist2Neighbors())

#In [436]: l = MergeList(l, v6.ObtainDist2Neighbors())

#In [437]: l = MergeList(l, v7.ObtainDist2Neighbors())

#In [438]: len(l)
#Out[438]: 606

#In [439]: i = bisect.bisect_left(l, v0)

#In [440]: i0 = bisect.bisect_left(l, v0)

#In [441]: i1 = bisect.bisect_left(l, v1)

#In [442]: i2 = bisect.bisect_left(l, v2)

#In [443]: i3 = bisect.bisect_left(l, v3)

#In [444]: i4 = bisect.bisect_left(l, v4)

#In [445]: i5 = bisect.bisect_left(l, v5)

#In [446]: i6 = bisect.bisect_left(l, v6)

#In [447]: i7 = bisect.bisect_left(l, v7)

#In [448]: c = Configuration(l)

#In [449]: c.LaplaceConstraint = np.zeros(len(l))

#In [450]: c.LaplaceConstraint[i0] = 1

#In [451]: c.LaplaceConstraint[i1] = -1

#In [452]: c.LaplaceConstraint[i2] = -1

#In [453]: c.LaplaceConstraint[i3] = 1

#In [454]: c.LaplaceConstraint[i4] = -1

#In [455]: c.LaplaceConstraint[i5] = 1

#In [456]: c.LaplaceConstraint[i6] = 1

#In [457]: c.LaplaceConstraint[i7] = -1

#In [458]: c.ObtainOptimizationValue()
#Out[458]: 0.014920089177509765

#In [459]: 0.014920089177509765/8
#Out[459]: 0.0018650111471887206

# This exceeds alpha

###########################################
# The second estimate starts here #########
###########################################

#########################################################################
# This case considers a node with distance 1 to P1, 1 to P2 and 2 to P3 #
#########################################################################

#In [460]: v0 = Vertex(0,0,0,0)

#In [461]: v1 = Vertex(1,1,0,0)

#In [462]: v2 = Vertex(2, -2, 0, 0)

#In [463]: v3 = Vertex(3, -1, 0, 0)

#In [464]: v4 = Vertex(0,0,1,1)

#In [465]: v5 = Vertex(1,1,1,1)

#In [466]: v6 = Vertex(2,-2,1,1)

#In [467]: v7 = Vertex(3,-1,1,1)

#In [468]: l = v0.ObtainDist2Neighbors()

#In [469]: l = MergeList(l, v1.ObtainDist3Neighbors())

#In [470]: l = v0.ObtainDist2Neighbors()

#In [471]: l = MergeList(l, v2.ObtainDist2Neighbors())

#In [472]: l = MergeList(l, v3.ObtainDist2Neighbors())

#In [473]: l = MergeList(l, v4.ObtainDist2Neighbors())

#In [474]: l = MergeList(l, v5.ObtainDist2Neighbors())

#In [475]: l = MergeList(l, v6.ObtainDist2Neighbors())

#In [476]: l = MergeList(l, v7.ObtainDist2Neighbors())

#In [477]: i0 = bisect.bisect_left(l,v0)

#In [478]: i1 = bisect.bisect_left(l,v1)

#In [479]: i2 = bisect.bisect_left(l,v2)

#In [480]: i3 = bisect.bisect_left(l,v3)

#In [481]: i4 = bisect.bisect_left(l,v4)

#In [482]: i5 = bisect.bisect_left(l,v5)

#In [483]: i6 = bisect.bisect_left(l,v6)

#In [484]: i7 = bisect.bisect_left(l,v7)

#In [485]: c = Configuration(l)

#In [486]: c.LaplaceConstraint = np.zeros(len(l))

#In [487]: c.LaplaceConstraint[i0]=1

#In [488]: c.LaplaceConstraint[i1]=-1

#In [489]: c.LaplaceConstraint[i2]=-1

#In [490]: c.LaplaceConstraint[i3]=1

#In [491]: c.LaplaceConstraint[i4]=-1

#In [492]: c.LaplaceConstraint[i5]=1

#In [493]: c.LaplaceConstraint[i6]=1

#In [494]: c.LaplaceConstraint[i7]=-1

#In [495]: c.ObtainOptimizationValue()
#Out[495]: 0.015487964150616554

#In [496]: 0.01548796415/8
#Out[496]: 0.00193599551875

# This exceeds alpha

##################################
# The third estimate starts here #
##################################

#####################################################################
# The third case has a node with distance 1 to P1, 2 to P2, 2 to P3 #
#####################################################################

#In [497]: v0
#Out[497]: (0,0,0,0) 

#In [498]: v1
#Out[498]: (1,1,0,0) 

#In [499]: v2
#Out[499]: (2,-2,0,0) 

#In [500]: v3
#Out[500]: (3,-1,0,0) 

#In [501]: v4 = Vertex(0,0,2,2)

#In [502]: v5 = Vertex(1,1,2,2)

#In [503]: v6 = Vertex(2,-2,2,2)

#In [504]: v7 = Vertex(3,-1,2,2)

#In [505]: l = v0.ObtainDist2Neighbors()

#In [506]: l = MergeList(l, v1.ObtainDist2Neighbors())

#In [507]: l = MergeList(l, v2.ObtainDist2Neighbors())

#In [508]: l = MergeList(l, v3.ObtainDist2Neighbors())

#In [509]: l = MergeList(l, v4.ObtainDist2Neighbors())

#In [510]: l = MergeList(l, v5.ObtainDist2Neighbors())

#In [511]: l = MergeList(l, v6.ObtainDist2Neighbors())

#In [512]: l = MergeList(l, v7.ObtainDist2Neighbors())

#In [513]: i0 = bisect.bisect_left(l, v0)

#In [514]: i1 = bisect.bisect_left(l, v1)

#In [515]: i2 = bisect.bisect_left(l, v2)

#In [516]: i3 = bisect.bisect_left(l, v3)

#In [517]: i4 = bisect.bisect_left(l, v4)

#In [518]: i5 = bisect.bisect_left(l, v5)

#In [519]: i6 = bisect.bisect_left(l, v6)

#In [520]: i7 = bisect.bisect_left(l, v7)

#In [521]: c = Configuration(l)

#In [522]: c.LaplaceConstraint = np.zeros(len(l))

#In [523]: c.LaplaceConstraint[i0] = 1

#In [524]: c.LaplaceConstraint[i1] = -1

#In [525]: c.LaplaceConstraint[i2] = -1

#In [526]: c.LaplaceConstraint[i3] = 1

#In [527]: c.LaplaceConstraint[i4] = -1

#In [528]: c.LaplaceConstraint[i5] = 1

#In [529]: c.LaplaceConstraint[i6] = 1

#In [530]: c.LaplaceConstraint[i7] = -1

#In [531]: c.ObtainOptimizationValue()
#Out[531]: 0.01659891756220965

#In [532]: 0.01659891756220965/8
#     ...: 
#     ...: 
#Out[532]: 0.0020748646952762064

# This exceeds alpha

##############################################################################################
# This concludes the check that the point appears at distance at most 2 from each hyperplane #
##############################################################################################


