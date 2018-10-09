#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 20:35:01 2018

@author: paulolomeu
"""

import numpy as np


def kernel_matrix_full(Nx, Nz, Npml, Cpml, omega, delta, rho, v):
    
    # propriedades indiretas
    b = 1./(rho)
    K = rho*(v**2)
    
    
    imag = 1j
    
    
    csix = np.zeros((Nx,1), dtype=np.complex128)

    for i in range(0,Nx):

        if(Npml > 1  and  i < Npml ): # dentro da primeira camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (i-1)/(Npml-1))

        elif(Npml > 1  and  i > (Nx-Npml)): # dentro da segunda camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (Nx-i)/(Npml-1))

        else: # dentro do modelo
          gammax = 0.0


        csix[i] = 1.0  + imag * gammax/omega



    csiz = np.zeros((Nx,1), dtype=np.complex128)

    for k in range(0,Nz):
        
        # so tem um else se so tiver a camada inferior
        if (Npml > 1  and  k > (Nz-Npml)): # dentro da segunda camada de absorcao horizontal (leito)
            gammaz = Cpml*np.cos(np.pi/2.0 * (Nz-k)/(Npml-1))

        elif (Npml > 1  and  k < Npml ): # dentro da primeira camada de absorcao vertical (topo)
            gammaz = Cpml*np.cos(np.pi/2.0 * (k-1)/(Npml-1))

        else: # dentro do modelo
            gammaz = 0.0


        csiz[k] = 1.0  + imag * gammaz/omega
        
    #print csiz
    
    
    #colptr = np.zeros((Nz*Nx,1))
    #rowind = np.zeros((Nz*Nx,1))
    #values = np.empty((Nz*Nx,1),dtype=np.complex128)
    
    
#    colptr = []
#    rowind = []
#    values = []
    
    
    #C = np.zeros((Nx*Nz,Nx*Nz),dtype=np.complex128)
    C = np.zeros((Nx*Nz,2*Nx*Nz))
    C = C.view(dtype=np.complex128)
    
#    print(C)
    
    
    
    
    for i in range(0,Nx):
    
        for j in range(0,Nz):
       
            linha = (i)*Nz + j   #numero do noh ij
            
            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
            col1 = linha	  # corresponde ao noh i,k    #numero do noh ij
            col2 = linha - 1  # corresponde ao noh i-1,k  #numero do noh da esquerda
            col3 = linha + 1  # corresponde ao noh i+1,k  #numero do noh da direita
            col4 = linha - Nz  # corresponde ao noh i,k-1  #numero do noh de cima
            col5 = linha + Nz  # corresponde ao noh i,k+1  #numero do noh de baixo
            
    #        print 'linha', linha
    #        print 'col1', col1
    #        print 'col2', col2
    #        print 'col3', col3
    #        print 'col4', col4
    #        print 'col5', col5
            
            
            # no interno - nao tem c.c. aplicada
            if(i != 0 and i != Nx-1 and j != 0 and j != Nz-1 ):
                
                # staggered grid five point stencil (hustedt et al., 2004 equacao 8)
                c2 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i-1,j] )/(csiz[j]+csiz[j-1])
                c3 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i+1,j] )/(csiz[j]+csiz[j+1])
                c4 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i,j-1] )/(csix[i]+csix[i-1])
                c5 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i,j+1] )/(csix[i]+csix[i+1])
                c1 = omega**2/K[i,j] - c2 -c3 - c4 - c5
                               
                
                C[linha,col1:col1+1] = c1
                C[linha,col2:col2+1] = c2
                C[linha,col3:col3+1] = c3
                C[linha,col4:col4+1] = c4
                C[linha,col5:col5+1] = c5            
                
            # borda esquerda 
            elif(i == 0):
                
                # confirmar se delta varia em x e z???????????????????????
#                c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
#                c1 = omega**2/K[i,j] -c3 
#                
#                C[linha,col1:col1+1] = c1 + LB
#                C[linha,col3:col3+1] = c3
                
                c1 = -(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
                c5 = +(1.0/delta)
                
                LB = c1
                LBI = c5
                
                C[linha,col1:col1+1] = LB
                C[linha,col5:col5+1] = LBI
    
                
            # borda direita    
            elif(i == Nx-1):
                
                # confirmar se delta varia em x e z???????????????????????
#                c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
#                c1 = omega**2/K[i,j] - c2         
#                
#                C[linha,col1:col1+1] = c1 + RB
#                C[linha,col2:col2+1] = c2
                
                c1 = +(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
                c4 = -(1.0/delta)
                
                RB = c1
                RBI = c4
                
                C[linha,col1:col1+1] = RB
                C[linha,col4:col4+1] = RBI   
                
                
            # topo
            elif(j == 0):
                
                # confirmar se delta varia em x e z???????????????????????
#                c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
#                c1 = omega**2/K[i,j] - c5            
#                
#                C[linha,col1:col1+1] = c1 + FS
#                C[linha,col5:col5+1] = c5    
                
                c1 = -(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
                c3 = +(1.0/delta)
                
                FS = c1
                FSI = c3
                
                C[linha,col1:col1+1] = FS
                C[linha,col3:col3+1] = FSI
                
            
            # fundo
            elif(j == Nz-1):
                
                # confirmar se delta varia em x e z???????????????????????
#                c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
#                c1 = omega**2/K[i,j] - c4       
#                
#                C[linha,col1:col1+1] = c1 + BB
#                C[linha,col4:col4+1] = c4
                
                c1 = +(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
                c2 = -(1.0/delta)
                
                BB = c1
                BBI = c2
                
                C[linha,col1:col1+1] = BB
                C[linha,col2:col2+1] = BBI
    
    
#    print(C)
    
    return C
    
    
    
    
