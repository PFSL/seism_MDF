#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 20:18:32 2018

@author: paulolomeu
"""

import numpy as np


def input_data_new_1():
        
    ## entradas p o programa
    # set x limits
    xmin = 0.
    xmax = 1000.
    
    # set z limits
    zmin = 0.
    zmax = 1000.
    
    # set number of elements in each axis
    Nx = 101
    Nz = 101
    
    # set deep of source and hydrophones
    deep = 500.
    
    # numero de hidrofones do sistema de aquisicao
    num_rec = 98
    
    # frequencias em Hz
    freq = [8]
    # amplitudes
    AA = [10]
    
    # propriedades de amortecimento
    Npml = 20. # numero de camadas de amortecimento
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
    for i in range(0,Nx):
        for j in range(0,Nz):
            grade[cont,0] = cont
            grade[cont,1] = x[i]
            grade[cont,2] = z[j]
            cont = cont + 1
            
    
    # physical properties
    # velocity 
    v = np.ones((Nx,Nz)) # matriz de velocidades
    v1 = 1500.
    v2 = 1500.
    v3 = 1500.
    
    # desity
    rho = np.ones((Nx,Nz)) # matriz de densidades (para meios heterogeneos)
    rho1 = 1000
    rho2 = 1000
    rho3 = 1000
    
#    rho_stk = np.array([rho1, rho2, rho3])
#    v_stk = np.array([v1, v2, v3])
    
    # matriz de velocidades de 3 camadas planas
    for i in range(0,Nx):
        for j in range(0,int(Nz)):
            rho[i,j] = rho[i,j]*rho1
            v[i,j] = v[i,j]*v1
        
#    for i in range(0,Nx):
#        for j in range(int(Nz/3), 2*int(Nz/3)):
#            rho[i,j] = rho[i,j]*rho2
#            v[i,j] = v[i,j]*v2
#    
#    for i in range(0,Nx):
#        for j in range(2*int(Nz/3),int(Nz)):
#            rho[j,i] = rho[j,i]*rho3
#            v[j,i] = v[j,i]*v3
    
    

#    print('VERIFICAÇAO DOS INPUTS DO PROBLEMA\n')
#    print('xmin, xmax - ',xmin, xmax)
#    print('zmin, zmax - ',zmin, zmax)
#    print('Nx, Nz     - ', Nx, Nz)
#    print('delta', delta,'\n')
#    print('profundidade da fonte - ', deep)
#    print('numero de receptores  - ', num_rec,'\n')
#    print('frequencias e amplitudes:')
#    for line in solicita:
#        print(line[0], line[1])
    
    
    return Nx, Nz, delta, grade, v, rho, Npml, Cpml, solicita, deep, num_rec, num_freq
    


