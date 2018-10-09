#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 23:47:03 2017

@author: paulolomeu
"""

import numpy as np


def input_data():
    
    
    ## entradas p o programa
    # set x limits
    xmin = 0.
    xmax = 1000.
    #print(xmin, xmax)
    
    # set z limits
    zmin = 0.
    zmax = 1000.
    #print(zmin, zmax)
    
    # set number of elements in each axis
    Nx = 3
    Nz = 3
    delta = (zmax - zmin)/(Nz-1)
    #print(Nx, Nz, delta)
    
    # create the grid
    x = np.linspace(xmin,xmax,Nx)
    z = np.linspace(zmin,zmax,Nz)
    #print(x)
    #print(z)
    # aloca espaco para a grade
    grade = np.zeros((Nx*Nz,3))
    #print(grade)
    
    cont = 0
    # monta lista com informacoes da grade
    for i in range(0,Nx):
        for j in range(0,Nz):
            grade[cont,0] = cont
            grade[cont,1] = x[i]
            grade[cont,2] = z[j]
            cont = cont + 1
    #print('grade')
            print(cont, x[i],z[j])
    
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
    
    rho_stk = np.array([rho1, rho2, rho3])
    v_stk = np.array([v1, v2, v3])
    
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
    
    #print(np.transpose(v)) # transposta para visualizaçao ficar conforme profundidade
    #print(np.transpose(rho)) # transposta para visualizaçao ficar conforme profundidade
     
    
    # propriedades indiretas
    b = 1./(rho)
    K = rho*(v**2)
    #print(np.transpose(b))
    #print(np.transpose(K))
    
    # propriedades de amortecimento
    Npml = 1. # numero de camadas de amortecimento
    Cpml = 50. # parametro obtido manualmente
    
    # elementos da solicitacao do sistema (tipo single-spread) *inserir split-spread depois
    #         * o o o o o o o o o o
    # o --> hidrofone   \\\   * --> air gun
    
    # numero de hidrofones do sistema de aquisicao
    num_hid = 10
    
    # frequencia em Hz
    freq = [24]
    # amplitude
    AA = [10] # amplitude do pulso
    
    
    #numero de tiros depende do tamanho do dominio
    #  CRIAR SISTEMA EM QUE A MATRIZ SOLICITA CONTEM TODOS OS TIROS COM POSIÇAO, FREQUENCIA E AMPLITUDE
    
    
    # vetor com informacoes das fontes solicitantes teste -> fonte no centro do dominio
    # [coord x, coord z, freq(Hz), ts, amplitude]
    solicita = [[500,500,freq[0],AA[0]]]
    
    
    solicita = np.array(solicita)
    #print(solicita)
    
    
    return Nx, Nz, delta, grade, v, rho, v_stk, rho_stk, b, K, Npml, Cpml, solicita

    


