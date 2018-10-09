#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 21:08:45 2017

@author: paulolomeu
"""

from MDF_kernel import k_matrix_sp_t

from solicita_vetor import solicita_vetor



def helmholtz_2d(Nx, Nz, grade, b, K, Npml, Cpml, omega, delta, v, solicita):
    
    
    Csp = k_matrix_sp_t(Nx,Nz,Npml,Cpml,omega,b,K,delta,v)
    
    
#    from solicita_vetor import *
    # montagem do vetor de solicitacao
    s, omega = solicita_vetor(Nx,Nz,solicita,grade)
    
    
    import scipy as sp
    psp = sp.sparse.linalg.spsolve(Csp,s)
    
    
    
    return psp