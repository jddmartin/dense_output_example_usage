"""
Integration of complex Lorenz equations to demonstrate dense output capabilities
of the DOP853 ode solver code.

Fowler, Gibbon and McGuinness,  "The complex Lorenz equations", 
Physica 4D, 139-163 (1982)
http://dx.doi.org/10.1016/0167-2789(82)90057-4
"""

import math, numpy
from scipy.integrate import ode, complex_ode, dense_dop
from collections import namedtuple

class ComplexLorenz(object):
    def __init__(self, b, sigma, r1, r2, e):
        self.b=b
        self.sigma=sigma
        self.r=r1+1.0j*r2
        self.a=1.0-1.0j*e

    def f(self, t, v):
        """Complex Lorenz system odes, Eq. 1.1 of Fowler et al. (see above).
        """
        (x, y, z) = v 
        return [-self.sigma*x+self.sigma*y,
                -x*z+self.r*x-self.a*y,
                -self.b*z+(x.conjugate()*y).real]

class SolOut(object):
    def __init__(self, tbounds=None):
        if tbounds is not None:
            Bounds=namedtuple("bounds", "min max")
            self.tbounds=Bounds(tbounds[0], tbounds[1])
        else:
            self.tbounds=None

    def solout(self, nr, told, t, v, con_view, icomp):
        x, y, z = v
        if nr > 1:
            if self.tbounds:
#                print self.tbounds.min, self.tbounds.max
                if ((t < self.tbounds.min) or (t > self.tbounds.max)):
                    return
            print t, y.real, x.real, z.real

# Figure 1 and 2:
#system=ComplexLorenz(4.0/3.0, 2.0, 2.0, 1.0, 3.0)
#aSolOut=SolOut()

# Figure 4:
system=ComplexLorenz(0.8, 2.0, 256.0, 1.0, 3.0)
aSolOut=SolOut(tbounds=(835.0, 935.0))

# initial conditions, and length of time to integrate:
x0=0.0
y0=[1.0, 0.0, 0.0]
xend=1000

# desired tolerances:
itol=0
rtol=1.0e-7
atol=rtol

system.f(0.0, [1, 2, 3])
ig = complex_ode(system.f).set_integrator('dopri5', atol=atol, rtol=rtol, 
                                          nsteps=1000000)

ig.set_solout(aSolOut.solout, dense_components=(0,1,))

rpar=[]
ig.set_initial_value(y0, x0) #.set_f_params(rpar)

ret = ig.integrate(xend)

iw=ig._integrator.returned_iwork
#print "     tol={0:5.2e}   fcn= {1:d} step= {2:d} accpt= {3:d} rejct= {4:d}".format(atol,*iw[16:20])
