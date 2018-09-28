#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 09:00:21 2017

@author:
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

# system variables
u = symbols('u')
U = [u]
parameters = []

#
# f_0(U)_t + F_1(U)_x = (B(U)U_x)_x
#

f0 = [u]

f1 = [u**2/2-u]

BU = [[1]]

G = [0]

# file path to location where code will be written
file_path = sys.path[0]
create_code(U,parameters,f0,f1,G,BU,file_path)

