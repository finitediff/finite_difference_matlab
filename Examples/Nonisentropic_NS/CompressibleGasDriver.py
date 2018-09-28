"""
Example of 1-dimensional Compressible Gas using an ideal gas law

Created on Thursday, August 3, 2017
@authors:  Blake Barker, Jalen Morgan, Taylor Paskett

see README.txt for help on using PyPDE.

************************************************************
*********************************************************"""

import pypde
import pde

if __name__ == "__main__":
    
    #p is initialized as a dictionary with default values.
    p = pypde.init()
    
    #Define the PDE
    p.update({
        #Define the scope of the PDE
        "T": 8,
        "L": 50,
        "pdeTPoints": 30,
        "pdeXPoints": 200,
        "fileName":"compressible_gas",
        
        #Define the PDE
        "pde": ["V_t + V_x - U_x",
                "U_t + U_x + p_x - (mu*U_x/V)_x",
                "E_t + E_x + (p*U)_x - (mu*U*U_x/V)_x - (k*T_x/V)_x "],
                
        "pdeUnknowns": ['V','U','e'],
        
        "pdeSubstitute": ["T = e/cnu", 
                          "p = Gamma*e/V", 
                          "E = e + U**2/2"],
        
        #Set the parameters for the PDE
        "cnu": 1.0,
        "Gamma": .666,
        "mu": 1.0,
        "k": 1.0,
        
        #Set files from which to load initial values
        "pdeInitialValueFiles": ["vProfile.txt", "uProfile.txt", "eProfile.txt"]
    })
    
    #Solve the system
    pde.solve(p)