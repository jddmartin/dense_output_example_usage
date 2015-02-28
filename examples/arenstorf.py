"""Calculation of an Arenstorf orbit illustrating usage of the
 dense output option of the DOPRI5 integrator.

An introductory discussion is given on pages 129-131 of Hairer et al.'s
"Solving Ordinary Differential Equations, Nonstiff Problems",
Second Revised Edition, Springer, 1993.
This code is modelled after "Driver for the code DORPI5" in the Appendix of
this book, available at: http://www.unige.ch/~hairer/prog/nonstiff/dr_dopri5.f
"""

from __future__ import division, print_function, absolute_import

import math
import numpy
from scipy.integrate import ode, complex_ode, dense_dop


def solout(nr, xold, x, y, con_view, icomp):
    global xout
    # Note that the use of "global" can be avoided using classes
    # (recommended). i.e. make solout a member function, and xout a class
    # attribute (see "complex_lorenz.py" as an example).

    # mimic Fortran output format:
    format_string = (" X = {0:5.2f}   Y ={1:18.10e} {2:18.10e}    "
                     "NSTEP = {3:4d}")
    if nr == 1:
        print(format_string.format(x, y[0], y[1], nr-1))
        xout += 2.0
    else:
        while x >= xout:
            dense = dense_dop(xout, xold, x, con_view)
            print(format_string.format(xout, dense[0], dense[1], nr-1))
            xout += 2.0


def f_arenstorf(x, y, rpar):
    """The system of differential equations.
    """
    amu, amup = rpar
    r1 = (y[0]+amu)**2+y[1]**2
    r1 = r1*math.sqrt(r1)
    r2 = (y[0]-amup)**2+y[1]**2
    r2 = r2*math.sqrt(r2)
    f2 = y[0]+2*y[3]-amup*(y[0]+amu)/r1-amu*(y[0]-amup)/r2
    f3 = y[1]-2*y[2]-amup*y[1]/r1-amu*y[1]/r2
    return [y[2], y[3], f2, f3]


# parameters for differential equation system:
rpar = numpy.zeros((2,), float)
rpar[0] = 0.012277471
rpar[1] = 1.0-rpar[0]

# initial conditions, and length of time to integrate:
x0 = 0.0
y0 = [0.994, 0.0, 0.0, -2.00158510637908252240537862224]
xend = 17.0652165601579625588917206249

# desired tolerances:
itol = 0
rtol = 1.0e-7
atol = rtol

xout = 0.0  # used to keep track of current time in "solout".

ig = ode(f_arenstorf).set_integrator('dopri5', atol=atol, rtol=rtol)

# Now use the new dense output option by specification of "dense_components"
# in call of ".set_solout" method.
# Although there are 4 components (2 positions and 2 velocities) we
# only request dense output for the 2 positions:
ig.set_solout(solout, dense_components=(0, 1,))

ig.set_initial_value(y0, x0).set_f_params(rpar)
ret = ig.integrate(xend)

print(" X = XEND    Y ={0:18.10e} {1:18.10e}".format(ret[0], ret[1],))
iw = ig._integrator.returned_iwork
print("     tol={0:5.2e}   fcn= {1:d} step= {2:d} "
      "accpt= {3:d} rejct= {4:d}".format(atol, *iw[16:20]))
