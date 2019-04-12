#/usr/bin/env python3
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
mu, nu  = symbols('mu, nu')

# system variables
u, v = symbols('u, v')

# vector containing the system variables
U = [u, v]

# system parameters
parameters = [mu, nu]

# flux form
#
# f_0(U)_t + F_1(U)_x = (B(U)U_x)_x
#

f0 = [u,v]

f1 = [0,0]

BU = [[0,-1],[1,0]]

G = [-mu*v + (u**2+v**2)*v, u-(u**2+v**2)*u+2*nu*v]


# file path to location where code will be written
file_path = sys.path[0]

create_code(U,parameters,f0,f1,G,BU,file_path)

