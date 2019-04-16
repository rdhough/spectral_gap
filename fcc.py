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
    def spherical_integrand(r,thet1, thet2):
        if r == 0:
            return 0
        else:
            return integrand(r*np.cos(thet1), r*np.sin(thet1)*np.cos(thet2),r*np.sin(thet1)*np.sin(thet2) )*r**2 * np.sin(thet1)
        
    return integrate.nquad(lambda r, thet1, thet2: spherical_integrand(r,thet1, thet2), [[0, .5],[0,np.pi], [0, 2*np.pi]])[0]

# Integration over the torus (R/Z)^3 with singularity at 0

def SingularityIntegral(integrand):
    def smooth_integrand(x,y,z):
        if x**2 + y**2+z**2 == 0:
            return 0
        else:
            return integrand(x,y,z)* (1-BumpFunction(4*(x**2 + y**2+z**2)))
    def singular_integrand(x,y,z):
        return integrand(x,y,z) * BumpFunction(4*(x**2 + y**2+ z**2))
    SmoothIntegral = integrate.nquad(lambda x,y,z: smooth_integrand(x,y,z), [[-.5,.5],[-.5,.5], [-.5,.5]])[0]
    SingularIntegral = SphericalCoordinateIntegrate(singular_integrand)
    return SmoothIntegral + SingularIntegral


def CosDif(x):
    return 1 - np.cos(2 *np.pi * x)


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

# This implements a point in the FCC lattice
class Vertex(object):
    def __init__(self, x_coor, y_coor, z_coor):
        self.x_coor = x_coor
        self.y_coor = y_coor
        self.z_coor = z_coor
    def __lt__(self, other):
        if self.x_coor < other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor < other.y_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor < other.z_coor:
            return True
        else:
            return False
    def __le__(self, other):
        if self.x_coor < other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor < other.y_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor <= other.z_coor:
            return True
        else:
            return False
    def __gt__(self, other):
        if self.x_coor > other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor > other.y_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor > other.z_coor:
            return True
        else:
            return False
    def __ge__(self, other):
        if self.x_coor > other.x_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor > other.y_coor:
            return True
        elif self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor >= other.z_coor:
            return True
        else:
            return False
    def __eq__(self, other):
        if self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor == other.z_coor:
            return True
        else:
            return False
    def __ne__(self, other):
        if self.x_coor == other.x_coor and self.y_coor == other.y_coor and self.z_coor == other.z_coor:
            return False
        else:
            return True
    def __getitem__(self, key):
        if key == 1:
            return self.x_coor
        if key == 2:
            return self.y_coor
        if key == 3:
            return self.z_coor
# Definition of 12 neighbors
    def ObtainNeighbors(self):
        v1 = Vertex(self.x_coor + 1, self.y_coor, self.z_coor)
        v2 = Vertex(self.x_coor - 1, self.y_coor, self.z_coor)
        v3 = Vertex(self.x_coor, self.y_coor+1, self.z_coor)
        v4 = Vertex(self.x_coor, self.y_coor-1, self.z_coor)
	v5 = Vertex(self.x_coor, self.y_coor, self.z_coor-1)
	v6 = Vertex(self.x_coor, self.y_coor, self.z_coor+1)
        v7 = Vertex(self.x_coor+1, self.y_coor-1, self.z_coor)
        v8 = Vertex(self.x_coor-1, self.y_coor+1, self.z_coor)
        v9 = Vertex(self.x_coor+1, self.y_coor, self.z_coor-1)
        v10 = Vertex(self.x_coor-1, self.y_coor, self.z_coor+1)
        v11 = Vertex(self.x_coor, self.y_coor-1, self.z_coor+1)
        v12 = Vertex(self.x_coor, self.y_coor+1, self.z_coor-1)
        return SortList([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12])
    
#Obtain the distance 2 neighborhood of a point
    def ObtainDist2Neighbors(self):
	l = self.ObtainNeighbors()
        returnlist = l
	for e in l:
	    returnlist = MergeList(returnlist, e.ObtainNeighbors())
        return returnlist

    def ObtainDist4Neighbors(self):
        l = self.ObtainDist2Neighbors()
        returnlist = l
        for e in l:
            returnlist = MergeList(returnlist, e.ObtainDist2Neighbors())
        return returnlist
#Obtain the distance 6 neighborhood of a point
    def ObtainDist6Neighbors(self):
        l = self.ObtainDist4Neighbors()
        returnlist = l
        for e in l:
            returnlist = MergeList(returnlist, e.ObtainDist2Neighbors())
        return returnlist
# Translate the point    
    def Transform(self, w1, w2, w3):
        x_coor = self.x_coor * (w1.x_coor) + self.y_coor * (w2.x_coor) + self.z_coor *(w3.x_coor)
        y_coor = self.x_coor * (w1.y_coor) + self.y_coor * (w2.y_coor) + self.z_coor *(w3.y_coor)
        z_coor = self.x_coor * (w1.z_coor) + self.y_coor * (w2.z_coor) + self.z_coor *(w3.z_coor)
        return Vertex(x_coor, y_coor, z_coor)
    def __str__(self):
        return "(" +str(self.x_coor) + "," +str(self.y_coor) +","+ str(self.z_coor)+ ") "
    def __repr__(self):
        return "(" +str(self.x_coor) + "," +str(self.y_coor) +","+ str(self.z_coor)+ ") "




    
class Configuration(object):
    def __init__(self, vertices):
        self.vertices = SortList(vertices)
        self.n = len(vertices)
	self.LaplaceConstraint = np.zeros(self.n)
	self.Laplacian = np.zeros(self.n)
        self.neighbors = []
        self.dist2neighbors = []
        for v in vertices:
            for n in v.ObtainNeighbors():
                if NotPresent(self.neighbors, n) and NotPresent(self.vertices, n):
                    bisect.insort(self.neighbors, n)
        for v in vertices:
            for n in v.ObtainDist2Neighbors():
                if NotPresent(self.dist2neighbors, n) and NotPresent(self.vertices, n):
                    bisect.insort(self.dist2neighbors, n)
        self.x_max = vertices[0].x_coor
        self.x_min = vertices[0].x_coor
        self.y_max = vertices[0].y_coor
        self.y_min = vertices[0].y_coor
        self.z_max = vertices[0].z_coor
        self.z_min = vertices[0].z_coor
        for ve in vertices:
            self.x_max = max(self.x_max, ve.x_coor)
            self.x_min = min(self.x_min, ve.x_coor)
            self.y_max = max(self.y_max, ve.y_coor)
            self.y_min = min(self.y_min, ve.y_coor)
            self.z_max = max(self.z_max, ve.z_coor)
            self.z_min = min(self.z_min, ve.z_coor)
        self.prevector = np.zeros((self.x_max - self.x_min + 1, self.y_max-self.y_min + 1, self.z_max - self.z_min+1))
        for ve in vertices:
            self.prevector[ve.x_coor-self.x_min, ve.y_coor-self.y_min, ve.z_coor - self.z_min] =1
        self.value = -1
        self.NumBins = 5
        v0 = Vertex(0,0,0)
        v0_neighbors = v0.ObtainNeighbors()
        self.rotation_list = []
        for w1 in v0_neighbors:
            joint_neighbors = []
            for v1 in w1.ObtainNeighbors():
                if not NotPresent(v0_neighbors, v1):
                    joint_neighbors.append(v1)
                for w2 in joint_neighbors:
                    for w3 in w2.ObtainNeighbors():
                        if not NotPresent(joint_neighbors, w3):
                            self.rotation_list.append([w1, w2, w3])
        
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

     

# Translate configurations by (x,y,z)

    def Translate(self, x, y, z):
        newverts = []
        for v in self.vertices:
            v_ = Vertex(v[1]+x, v[2]+y, v[3] + z)
            newverts.append(v_)
        return Configuration(SortList(newverts))


    def Transform(self, w1, w2, w3):
        newverts = []
        for v in self.vertices:
	    newverts.append(v.Transform(w1,w2,w3))
        return Configuration(SortList(newverts))

    def Normalize(self):
        cc = self.Translate(-self.x_min, -self.y_min, -self.z_min)
        return cc


    def NumVariables(self):
        return len(self.vertices) + len(self.neighbors)

    def NumVertexVariables(self):
        return len(self.vertices)

    def setLaplacian(self, v):
	self.Laplacian = v

    def setLaplaceConstraint(self, v):
	self.LaplaceConstraint = v

# Sum of squares function
        
    def ObtainFunctionValue(self, x):
        retval = 0
        for y in x:
            retval +=  y*y
        return retval

# Gradient for sum of squares    
    
    def ObtainDerivative(self, x):
        deriv = np.zeros(self.NumVariables())
        for i in range(self.NumVariables()):
            deriv[i] = 2 * x[i]
        return deriv

# Hessian for sum of squares
    
    def ObtainHessian(self, x):
        Hessian = np.zeros((self.NumVariables(), self.NumVariables()))
        for i in range(self.NumVariables()):
            Hessian[i,i] = 2 
        return Hessian
    
# Signed bounds for Q'
    def ObtainBounds(self):
        l=[]
        for i in range(self.NumVariables()):
            l.append((-.5, .5))
        return tuple(l)
   
# Signed constraint matrix for Q'
    def ObtainConstraintMatrix(self):
        nv = self.NumVariables()
        nvv = self.NumVertexVariables()
        returnmat = np.zeros((nvv, nv))
        for i in range(nvv):
            returnmat[i,i] = 12
            for n in self.vertices[i].ObtainNeighbors():
                if not NotPresent(self.vertices, n):
                    ind = bisect.bisect(self.vertices, n)-1
                else:
                    ind = nvv + bisect.bisect(self.neighbors, n)-1
                returnmat[i, ind]=-1
        return returnmat
            

# This method implements the optimization program Q'
    
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

            
    def CheckValue(self):
        if self.value < 0:
            return self.ObtainValue()
        else:
            return self.value

# The L2 norm calculated by Parseval    
    def ObtainL2(self):
        def integrand(x,y,z):
            rv = 0
            for k in range(self.n):
                rv += self.Laplacian[k] * cmath.exp(-2 *np.pi *1j *(self.vertices[k].x_coor*x + self.vertices[k].y_coor*y + self.vertices[k].z_coor*z))
            return abs(rv)**2 / (12 - 2.*(np.cos(2*np.pi *x)+ np.cos(2 * np.pi*y) + np.cos(2*np.pi*z) +np.cos(2*np.pi * (x-y))+np.cos(2*np.pi * (x-z)) + np.cos(2*np.pi * (y-z))))**2
        
        return  SingularityIntegral(integrand)

# The value of Xi is determined for n_1v_1 + n_2v_2 + n_3v_3 such that |n_i| <= n
    def ObtainXiRepn(self, n):
        def FT(x,y,z):
            rv = 0
            for k in range(self.n):
                rv += self.Laplacian[k] * cmath.exp(-2 *np.pi *1j *(self.vertices[k].x_coor*x + self.vertices[k].y_coor*y + self.vertices[k].y_coor*z))
            return rv / (12 - 2.*(np.cos(2*np.pi *x)+ np.cos(2 * np.pi*y) + np.cos(2*np.pi*z) +np.cos(2*np.pi * (x-y))+np.cos(2*np.pi * (x-z)) + np.cos(2*np.pi * (y-z))))
        time_domain = np.zeros((2*n+1, 2*n+1, 2*n+1))
        for a in range(2*n+1):
            for b in range(2*n+1):
                for c in range(2*n+1):
                    def integrand(x,y,z):
                        rvalue = FT(x,y,z) * cmath.exp(-2 * np.pi *1j *((n-a)*x + (n-b)*y + (n-c)*z))
                        return rvalue.real
                    time_domain[a,b,c] =  SingularityIntegral(integrand)
        return time_domain

# Sums the functional in a neighborhood of 0, then makes a correction using the L2 norm on the complement.
    
    def ObtainValue(self, n):
        self.value = 0
        l2 = self.ObtainL2()
        local_l2 = 0
        Xi = self.ObtainXiRepn(n)
        for a in range(2*n+1):
            for b in range(2*n+1):
                for c in range(2*n+1):
                    self.value += 1 - np.cos(2 * np.pi * Xi[a,b,c])
                    local_l2 += Xi[a,b,c]**2
        excess_l2 = l2 - local_l2
        self.value += 2*np.pi**2 * excess_l2 - (1./3.)*np.pi**4 *excess_l2**2
        error = (1./3.)*np.pi**4 *excess_l2**2
        return [self.value, error]

# In the plane, applying a rotation by 60 several times it is possible to obtain a position with non-negative v1, v2 coordinates, similarly v3.  Then apply a permutation to make v1 >=v2>=v3.
    
def SmallXiSearch():
    v0 = Vertex(0,0,0)
    l = v0.ObtainDist6Neighbors()
    test_list = []
    for v in l:
        if v != v0  and v[1]>= v[2] and v[2] >= v[3] and v[3] >=0:
            c = Configuration([v, v0])
            c.Laplacian = [1,-1]
            val = c.ObtainL2()
            print [v,val]
            if val < 0.01963:
                test_list.append(c)
    return test_list
            
    
#In [17]: c = Configuration([v0,v1])

#In [18]: c.Laplacian = [1,-1]

#In [19]: c.ObtainValue(5)
#Out[19]: [0.36239745177151528, 2.0519394474897139e-05]
   




#In [68]: c = Configuration([v])

#In [69]: c.LaplaceConstraint = [1]

#In [71]: c.ObtainOptimizationValue()
#Out[71]: 0.006410256422361247

#In [72]: l1
#Out[72]: 
#[(-1,0,0) ,
# (-1,0,1) ,
# (-1,1,0) ,
# (0,-1,0) ,
# (0,-1,1) ,
# (0,0,-1) ,
# (0,0,0) ,
# (0,0,1) ,
# (0,1,-1) ,
# (0,1,0) ,
# (1,-1,0) ,
# (1,0,-1) ,
# (1,0,0) ]
#In [73]: c1 = Configuration(l1)

#In [74]: c1.LaplaceConstraint = np.zeros(13)
#In [76]: c1.LaplaceConstraint[6]=1

#In [78]: c1.ObtainOptimizationValue()
#Out[78]: 0.009579728060137939

#In [79]: l2 = v.ObtainDist2Neighbors()

#In [80]: c2 = Configuration(l2)

#In [81]: c2.LaplaceConstraint = np.zeros(55)

#In [82]: c2.LaplaceConstraint[27] =1

#In [83]: c2.vertices[27]
#Out[83]: (0,0,0) 

#In [84]: c2.ObtainOptimizationValue()
#Out[84]: 0.012581342693839913

# This check verifies that the best configuration is two adjacent nodes with opposite signs

#In [117]: SmallXiSearch()
#[(1,0,0) , 0.01867584976614521]
#[(1,1,0) , 0.030287569856061142]
#[(1,1,1) , 0.04083921502222432]
#[(2,0,0) , 0.03422009034146848]
#[(2,1,0) , 0.04370401666881218]
#[(2,1,1) , 0.05334810022634776]
#[(2,2,0) , 0.055477360305846504]
#[(2,2,1) , 0.06486516148302635]
#[(2,2,2) , 0.0758894762889398]
#[(3,0,0) , 0.04880606985568302]
#[(3,1,0) , 0.05749595022007078]
#[(3,1,1) , 0.06657010458739787]
#[(3,2,0) , 0.06823424525968029]
#[(3,2,1) , 0.07732691326640999]
#[(3,3,0) , 0.08011841765111065]
#[(4,0,0) , 0.06312173880108815]
#[(4,1,0) , 0.07141102890551455]
#[(4,1,1) , 0.08011165213272609]
#[(4,2,0) , 0.0814691025383955]
#[(5,0,0) , 0.0773350384581015]
#[(5,1,0) , 0.08538143925970967]
#[(6,0,0) , 0.09149888929233442]

