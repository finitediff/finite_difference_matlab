#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 09:00:21 2017

@author: blake
"""

# -----------------------------------------------------------------------------
# imports
# -----------------------------------------------------------------------------

from sympy import symbols
from time_evolution import create_code
from sympy import exp

# -----------------------------------------------------------------------------
# setup
# -----------------------------------------------------------------------------

# parameters
nu, kappa, cnu, A, q, k, G, Ti, spd = symbols('nu, kappa, cnu, A, q, k, G, Ti, spd')

# system variables
tau, u, e, z = symbols('tau, u, e, z')

# vector containing the system variables
U = [tau, u, e, z]

# system parameters
parameters = [nu, kappa, D, cnu, A, q, k, G]

#
# f_0(U)_t + F_1(U)_x = (B(U)U_x)_x
#

f0 = [tau, u, e+u**2/2,z]

f1 = [-u-spd*tau, G*e/tau-spd*u, G*e*u/tau-spd*(e+u**2/2),-spd*z]

BU = [[0,0,0,0],[0,nu/tau,0,0],[0,nu*u/tau,kappa/(cnu*tau),0],[0,0,0,D/tau**2]]

G = [0,0,-q*k*exp(-A/(e/cnu-Ti))*z,k*exp(-A/(e/cnu-Ti))*z]


# file path to location where code will be written
file_path = '/Users/blake/Dropbox/stablab20/rNS/matlab'

create_code(U,parameters,f0,f1,G,BU,file_path)

