"""
Integration of complex Lorenz equations to demonstrate dense output 
capabilities of the DOP853 ode solver code.

The goal of this example is to reproduce Fig. 1 from
  Fowler, Gibbon and McGuinness,  "The complex Lorenz equations", 
  Physica 4D, 139-163 (1982)
  http://dx.doi.org/10.1016/0167-2789(82)90057-4
showing the advantage of dense output.
"""

from __future__ import division, print_function, absolute_import

import math, numpy
from scipy.integrate import ode, complex_ode, dense_dop
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt

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
    def __init__(self, tbounds=None, tinc=0.1):
        if tbounds is not None:
            Bounds=namedtuple("bounds", "min max")
            self.tbounds=Bounds(tbounds[0], tbounds[1])
        else:
            self.tbounds=None
        self.sparse_output=[]
        self.dense_output=[]
        self.tinc=tinc

    def solout(self, nr, told, t, v, con_view, icomp):
        x, y, z = v
        if nr == 1: # initial starting point:
            self.sparse_output.append((t, x.real, y.real, z.real,))
            self.dense_output.append((t, x.real, y.real, z.real,))
            self.tdense=t+self.tinc
        else: # subsequent step positions (after starting point):
            if self.tbounds:
                if ((t < self.tbounds.min) or (t > self.tbounds.max)):
                    return
            self.sparse_output.append((t,x.real,y.real,z.real,))
            while t > self.tdense:
                xd, yd, zd = dense_dop(self.tdense, told, t, con_view)
                self.dense_output.append((self.tdense,xd.real,yd.real,zd.real))
                self.tdense += self.tinc

if __name__ == "__main__":
    # Figure 1 and 2:
    system=ComplexLorenz(4.0/3.0, 2.0, 2.0, 1.0, 3.0)
    aSolOut=SolOut(tinc=0.01)

    # initial conditions, and length of time to integrate:
    t0=0.0
    v0=[1.0, 0.0, 0.0]
    tend=10.0

    # desired tolerances:
    itol=0
    rtol=1.0e-7
    atol=rtol

    ig = complex_ode(system.f).set_integrator('dop853', atol=atol, rtol=rtol, 
                                              nsteps=1000000)
    ig.set_solout(aSolOut.solout, dense_components=(0,1,2))
    ig.set_initial_value(v0, t0)

    # solve system:
    ret = ig.integrate(tend)

    # plot dense and sparse outputs, emulating Fowler et al.'s Fig. 1.
    dense_output_array=np.asarray(aSolOut.dense_output)
    sparse_output_array=np.asarray(aSolOut.sparse_output)
    for plotnum, second_dataset in ((1,sparse_output_array),
                                    (2,dense_output_array)):
        plt.subplot(2, 1, plotnum, aspect=3.6/2.1, 
                    xlim=[-1.8,1.8], ylim=[-0.9,1.20])
        plt.plot(sparse_output_array[:,2], sparse_output_array[:,1],'+')
        plt.plot(second_dataset[:,2], second_dataset[:,1])
        plt.ylabel("real(X)") # unconventional X and Y ordering in paper!
        plt.xticks(np.arange(-1,2, 1.0))
        plt.minorticks_on()
        for point in (0,-1):
            plt.text(sparse_output_array[point,2],
                     sparse_output_array[point,1],
                     ('t=%.2f' % sparse_output_array[point,0]),fontsize=8,
                     verticalalignment='bottom', horizontalalignment='center')
    plt.xlabel("real(Y)")

    plt.savefig('reproduction_of_figure_1_of_fowler_et_al.png')
    plt.show()
