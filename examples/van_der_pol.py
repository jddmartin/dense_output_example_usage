"""Solution of van der Pol differential equation illustrating usage of the
 dense output option of the DOP853 integrator.

An introductory discussion is given on pages 111 of Hairer et al.'s
"Solving Ordinary Differential Equations, Nonstiff Problems",
Second Revised Edition, Springer, 1993.
but epsilon has a different definition here (reciprocal of the book's usage).

This code is modelled after "Driver for DOP853 on Van der Pol's equation"
available at: http://www.unige.ch/~hairer/prog/nonstiff/dr_dop853.f
(note that Hairer's source code comments accidentally mislabels this
file as a driver for DOPRI5, rather than DOP853).

Since being incorporated into scipy, the DOP853 code has been updated by Hairer
to have an "iout=3" dense option, which eliminates wasteful calculation of
interpolation polynomials in regimes where dense output is not required.
So although the following code reproduces the output of "dr_dop853.py", it
does so with more function evaluations (compare the last line of output with
Fortran code that uses "iout=3").
"""

from __future__ import division, print_function, absolute_import

import math
import numpy
from scipy.integrate import ode, complex_ode, dense_dop


def solout(nr, xold, x, y, con):
    global xout
    # Note that the use of "global" can be avoided using classes
    # (recommended). i.e. make solout a member function, and xout a class
    # attribute (see "complex_lorenz.py" as an example).

    # mimic Fortran output format:
    format_string = (" X = {0:5.2f}   Y ={1:18.10e}{2:18.10e}    "
                     "NSTEP = {3:4d}")
    if nr == 1:
        print(format_string.format(x, y[0], y[1], nr-1))
        xout = 0.1
    else:
        while x >= xout:
            dense = dense_dop(xout, xold, x, con)
            print(format_string.format(xout, dense[0], dense[1], nr-1))
            xout += 0.1


def f_van_der_pol(x, y, rpar):
    """The system of differential equations.
    """
    eps = rpar[0]
    return [y[1], ((1-y[0]**2)*y[1]-y[0])/eps]

# parameters for differential equation system:
rpar = numpy.zeros((1,), float)
rpar[0] = 1.0e-3

# initial conditions, and length of time to integrate:
x0 = 0.0
y0 = [2.0, 0.0]
xend = 2.0

# desired tolerances:
itol = 0
rtol = 1.0e-9
atol = rtol

xout = 0.0  # used to keep track of current time in "solout".

ig = ode(f_van_der_pol).set_integrator('dop853', atol=atol, rtol=rtol,
                                       nsteps=10000)

# Now use the new dense output option by specification of "dense_components"
# in call of ".set_solout" method:
ig.set_solout(solout, dense_components=(0, 1,))

ig.set_initial_value(y0, x0).set_f_params(rpar)
ret = ig.integrate(xend)

print(" X = XEND    Y ={0:18.10e}{1:18.10e}".format(ret[0], ret[1],))

iw = ig._integrator.returned_iwork
print("        tol={0:5.2e} \n fcn= {1:d} step= {2:d} "
      "accpt= {3:d} rejct= {4:d}".format(atol, *iw[16:20]))
