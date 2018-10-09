#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 20:18:36 2018

@author: paulolomeu
"""

import numpy as np


def input_data_new_3():
        
    ## entradas p o programa
    # set x limits
    xmin = 0.
    xmax = 1000.
    
    # set z limits
    zmin = 0.
    zmax = 1000.
    
    # set number of elements in each axis
    Nx = 11
    Nz = 11
    
    # set deep of source and hydrophones
    deep = 75.
    
    # numero de hidrofones do sistema de aquisicao
    num_rec = 6
    
    # frequencias em Hz
    freq = [2, 4]
    # amplitudes
    AA = [10, 20]
    
    # propriedades de amortecimento
    Npml = 1. # numero de camadas de amortecimento
    Cpml = 50. # parametro obtido manualmente
    
    
    # numero de frequencias avaliadas
    num_freq = np.size(freq)
#    print('numero de freq', num_freq)
    # martiz de frequencias e amplitudes
    solicita = []
    
    for src in range(0, num_freq):
        solicita.append(np.array([freq[src],AA[src]]))
    
    
    delta = (zmax - zmin)/(Nz-1)
    
    # create the grid
    x = np.linspace(xmin,xmax,Nx)
    z = np.linspace(zmin,zmax,Nz)
    
    # aloca espaco para a grade
    grade = np.zeros((Nx*Nz,3))
    
    cont = 0
    # monta lista com informacoes da grade
    for i in range(0,Nz):
        for j in range(0,Nx):
            grade[cont,0] = cont
            grade[cont,1] = x[j]
            grade[cont,2] = z[i]
            cont = cont + 1
            
    
    # physical properties
    # velocity 
    v = np.ones((Nx,Nz)) # matriz de velocidades
    v1 = 1500.
    v2 = 1600.
    v3 = 1700.
    
    # desity
    rho = np.ones((Nx,Nz)) # matriz de densidades (para meios heterogeneos)
    rho1 = 1000
    rho2 = 1100
    rho3 = 1200
    
#    rho_stk = np.array([rho1, rho2, rho3])
#    v_stk = np.array([v1, v2, v3])
    
    # matriz de velocidades de 3 camadas planas
    for i in range(0,Nx):
        for j in range(0,int(Nz/3)):
            rho[j,i] = rho[j,i]*rho1
            v[j,i] = v[j,i]*v1
        
    for i in range(0,Nx):
        for j in range(int(Nz/3), 2*int(Nz/3)):
            rho[j,i] = rho[j,i]*rho2
            v[j,i] = v[j,i]*v2
    
    for i in range(0,Nx):
        for j in range(2*int(Nz/3),int(Nz)):
            rho[j,i] = rho[j,i]*rho3
            v[j,i] = v[j,i]*v3
    
    
    
    return Nx, Nz, delta, grade, v, rho, Npml, Cpml, solicita, deep, num_rec, num_freq
    

