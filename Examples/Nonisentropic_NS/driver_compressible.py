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
import sys, os
sys.path.insert(1,os.path.dirname(os.path.dirname(sys.path[0]))+'/Core')
from time_evolution import create_code
from sympy import exp
import os

# -----------------------------------------------------------------------------
# setup
# -----------------------------------------------------------------------------

# parameters
nu, mu, eta, cnu, G = symbols('nu, mu, eta, cnu, G')

# system variables
rho, u, e = symbols('rho, u, e')

# vector containing the system variables
U = [rho,u,e]

# system parameters
parameters = [nu, mu, eta, cnu, G]

# flux form
#
# f_0(U)_t + F_1(U)_x = (B(U)U_x)_x
#

f0 = [rho,rho*u,rho*(e+u**2/2)]

f1 = [rho*u, rho*u**2+G*rho*e, rho*u*(e+u**2/2)+G*rho*u*e]

BU = [[0,0,0],[0,(2*mu+eta),0],[0,(2*mu+eta)*u, nu]]

G = [0,0,0]


# file path to location where code will be written
file_path = sys.path[0]

create_code(U,parameters,f0,f1,G,BU,file_path)

