#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 21:56:59 2017

@author: blake
"""

# -----------------------------------------------------------------------------
# code
# -----------------------------------------------------------------------------

# f0(U)_t+f1(U)_x + G(U)= (B(U)*U_x)_x


#
# left multiply a vector by a matrix
#

from sympy import diff
from sympy import symbols
from sympy import sympify

def matrix_vector_multiply(A,v):
    
    out = []
    for j in range(0,len(A)):
        temp = sympify(0)
        for k in range(0,len(A[0])):
            temp = temp+A[j][k]*v[k]
        out.append(temp)
        
    return out

#
# Jacobian
# 

def compute_jacobian(F,U):
    
    DF = []
    for j in range(0,len(U)):
        temp = []
        for k in range(0,len(U)):
            temp.append(diff(F[j],U[k]))
        DF.append(temp)
        
    return DF

# write the parameters in the matlab file
def write_common(fid,parameters,U):
        
    for j in range(0,len(parameters)):
        fid.write(str(parameters[j])+' = p.'+str(parameters[j])+';\n')
    
    for j in range(0,len(U)):
            
        fid.write(str(U[j]))
        fid.write('_o = U_o('+str(j+1)+',2:end-1);\n')
        fid.write(str(U[j]))
        fid.write('_om = U_o('+str(j+1)+',1:end-2);\n')
        fid.write(str(U[j]))
        fid.write('_op = U_o('+str(j+1)+',3:end);\n')
        
        fid.write(str(U[j]))
        fid.write('_n = U_n('+str(j+1)+',2:end-1);\n')
        fid.write(str(U[j]))
        fid.write('_nm = U_n('+str(j+1)+',1:end-2);\n')
        fid.write(str(U[j]))
        fid.write('_np = U_n('+str(j+1)+',3:end);\n')

#
# Rest of code
#


def create_code(U,parameters,f0,f1,H,BU,file_path):

    # get the first and second derivatives U_x and U_xx and time derivative U_t
    U_x = []
    U_xx = []
    U_t = []
    for j in range(0,len(U)):
        U_x.append(symbols(str(U[j])+"_x"))
        U_xx.append(symbols(str(U[j])+"_xx"))
        U_t.append(symbols(str(U[j])+"_t"))
        
    
    
    # Jacobian of f0
    Df0 = compute_jacobian(f0,U)
    
    # Jacobian of f1
    Df1 = compute_jacobian(f1,U)
    
    # Jacobian of G
    #DG = compute_jacobian(G,U)
        
    # get the derivative of B(U)
    dBU = []
    for j in range(0,len(U)):
        row = []
        for k in range(0,len(U)):
            temp = sympify(0)
            for l in range(0,len(U)):
                temp = temp+diff(BU[j][k],U[l])*U_x[l]
            row.append(temp)
        dBU.append(row)
        
    
    df0_Ut = matrix_vector_multiply(Df0,U_t)
        
    df1_Ux = matrix_vector_multiply(Df1,U_x)
    
    BU_Uxx = matrix_vector_multiply(BU,U_xx)
    
    dBU_Ux = matrix_vector_multiply(dBU,U_x)
    
    # The function
    F = []
    for j in range(0,len(df0_Ut)):
        F.append(df0_Ut[j]+df1_Ux[j]+H[j]-BU_Uxx[j]-dBU_Ux[j])
    
    # substitute derivatives
    
    Wo = []
    Wo_p = []
    Wo_m = []
    Wn = []
    Wn_p = []
    Wn_m = []
    for j in range(0,len(U)):
        # U_j^n
        Wo.append(symbols(str(U[j])+'_o'))
        # U_{j+1}^n
        Wo_p.append(symbols(str(U[j])+'_op'))
        # U_{j-1}^n
        Wo_m.append(symbols(str(U[j])+'_om'))
        # U_{j}^{n+1}
        Wn.append(symbols(str(U[j])+'_n'))
        # U_{j+1}^{n+1}
        Wn_p.append(symbols(str(U[j])+'_np'))
        # U_{j-1}^{n+1}
        Wn_m.append(symbols(str(U[j])+'_nm'))
        
    K,H = symbols('K,H')
        
    fd_Ut = []
    for j in range(0,len(U)):
        fd_Ut.append((Wn[j]-Wo[j])/K)
    
    fd_Ux = []
    for j in range(0,len(U)):
        fd = (Wn_p[j]-Wn_m[j])/(4*H)+(Wo_p[j]-Wo_m[j])/(4*H)
        fd_Ux.append(fd)
        
    fd_Uxx = []
    for j in range(0,len(U)):
        fd = (Wn_p[j]-2*Wn[j]+Wn_m[j])/(2*H**2)+(Wo_p[j]-2*Wo[j]+Wo_m[j])/(2*H**2)
        fd_Uxx.append(fd)
     
    # F converted to finite difference
    G = []
    for j in range(0,len(U)):
        temp = F[j]
        for k in range(0,len(U)):
            temp = temp.subs(U[k],Wn[k])
            temp = temp.subs(U_t[k],fd_Ut[k])
            temp = temp.subs(U_x[k],fd_Ux[k])
            temp = temp.subs(U_xx[k],fd_Uxx[k])
        G.append(temp)
            
    # Jacobian of G, a 3 dimensional array with first component corresponding to the 
    # equation, second component corresponding to the system variable, and third 
    # component corresponding to the backward, forward, or current spatial node
    Jac_G = []
    for j in range(0,len(U)):
        U_J = []
        for k in range(0,len(U)):
            temp = []
            temp.append(diff(G[j],Wn_m[k]))
            temp.append(diff(G[j],Wn[k]))
            temp.append(diff(G[j],Wn_p[k]))
            U_J.append(temp)
        Jac_G.append(U_J)
        
    
    # open file for writing
    fid = open(file_path+'/'+'fd_jac.m', "w+")
    	
    # introductory part
    fid.write('function out = fd_jac(U_n,U_o,K,H,p)')
    fid.write('\n%\n%Jacobian of finite difference scheme\n\n')
    fid.write('\n\n%Jacobian\n\n')
    
    write_common(fid,parameters,U)
    
    # write the equations
    for j in range(0,len(Jac_G)):
        for k in range(0,len(Jac_G[j])):
            for l in range(0,len(Jac_G[j][k])):
                fid.write('%\nout{'+str(j+1)+'}{'+str(k+1)+'}{'+str(l+1)+'} = ')
                temp = str(Jac_G[j][k][l])
                temp = temp.replace('**','.^')
                temp = temp.replace('*','.*')
                temp = temp.replace('/','./')
                fid.write(temp)
                fid.write(';\n')
                
    fid.close()
    
    
    #
    # write the function file
    #
    
    # open file for writing
    fid = open(file_path+'/'+'fd_F.m', "w+")
    	
    # introductory part
    fid.write('function out = fd_F(U_n,U_o,K,H,p)')
    fid.write('\n%\n%Function of the finite difference scheme\n\n')
    
    write_common(fid,parameters,U)
    
    fid.write('\n\nout = [')
    
    for j in range(0,len(G)):
        
        temp = str(G[j])
        temp = temp.replace('**','.^')
        temp = temp.replace('*','.*')
        temp = temp.replace('/','./')
        fid.write(temp)
        fid.write(';\n')
    
    fid.write('];\n\n')
    
    fid.close()
