#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 20:26:20 2018

@author: paulolomeu
"""

import numpy as np


#from numba import jit
#
#@jit
def kernel_matrix_sparse_opt_jit(Nx, Nz, Npml, Cpml, omega, delta, rho, v):
    
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


    """
    CONCEITO DE ESPARSIDADE VAI SER REDUZIDO A UMA MATRIZ PENTA BANDA COM OS
    TERMOS REPETIDOS COM COORDENADAS IGUAIS A ZERO E VALOR IGUAL A ZERO, LOGO
    NA SOMA DELES NAO HAVERA PROBLEMA NA MONTAGEM DA MATRIZ CSR
    """
    from scipy.sparse import csr_matrix    
    row = []
    col = []
    data = []
    
    row  = np.zeros((5*Nx*Nz), dtype=np.int)
    col  = np.zeros((5*Nx*Nz), dtype=np.int)
    data = np.zeros((5*Nx*Nz,1), dtype=np.complex128)


    cont = 0
    
    # no interno - nao tem c.c. aplicada
    for i in range(1,Nx-1):

        for j in range(1,Nz-1):

            linha = (i)*Nz + j   #numero do noh ij
            
            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
            col1 = linha	  # corresponde ao noh i,k    #numero do noh ij
            col2 = linha - 1  # corresponde ao noh i-1,k  #numero do noh da esquerda
            col3 = linha + 1  # corresponde ao noh i+1,k  #numero do noh da direita
            col4 = linha - Nz  # corresponde ao noh i,k-1  #numero do noh de cima
            col5 = linha + Nz  # corresponde ao noh i,k+1  #numero do noh de baixo


            # staggered grid five point stencil (hustedt et al., 2004 equacao 8)
            c2 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i-1,j] )/(csiz[j]+csiz[j-1])
            c3 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i+1,j] )/(csiz[j]+csiz[j+1])
            c4 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i,j-1] )/(csix[i]+csix[i-1])
            c5 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i,j+1] )/(csix[i]+csix[i+1])
            c1 = omega**2/K[i,j] - c2 -c3 - c4 - c5
                           
            
            row[cont]  = linha
            col[cont]  = col1
            print('c1 ',c1[0])
            data[cont,0] = c1[0]
            
            cont = cont + 1
            
#            row.append(linha)
#            col.append(col1)
#            data.append(np.complex128(c1))
            
            row[cont]  = linha
            col[cont]  = col2
            data[cont] = c2
            
            cont = cont + 1
            
#            row.append(linha)
#            col.append(col2)
#            data.append(np.complex128(c2))
            
            row[cont]  = linha
            col[cont]  = col3
            data[cont] = c3
            
            cont = cont + 1
            
#            row.append(linha)
#            col.append(col3)
#            data.append(np.complex128(c3))
            
            row[cont]  = linha
            col[cont]  = col4
            data[cont] = c4
            
            cont = cont + 1
            
#            row.append(linha)
#            col.append(col4)
#            data.append(np.complex128(c4))
            
            row[cont]  = linha
            col[cont]  = col5
            data[cont] = c5
            
            cont = cont + 1
            
#            row.append(linha)
#            col.append(col5)
#            data.append(np.complex128(c5))
            
            
            
            
            
            
   
    
    for j in range(0,Nz):
        
        # borda esquerda
        i = 0
         
        linha = (i)*Nz + j   #numero do noh ij
        
        col1 = linha	  # corresponde ao noh i,k    #numero do noh ij
        col5 = linha + Nz  # corresponde ao noh i,k+1  #numero do noh de baixo

            
        c1 = -(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
        c5 = +(1.0/delta)
        
        LB = c1
        LBI = c5

#        row.append(linha)
#        col.append(col1)
#        data.append([np.complex128(LB)])
        
        row[cont]  = linha
        col[cont]  = col1
        data[cont] = LB
        
        cont = cont + 1
            
        
#        row.append(linha)
#        col.append(col5)
#        data.append([np.complex128(LBI)])

        row[cont]  = linha
        col[cont]  = col5
        data[cont] = LBI
        
        cont = cont + 1
            

        # borda direita
        i = Nx-1
        
        linha = (i)*Nz + j   #numero do noh ij
        
        col1 = linha	  # corresponde ao noh i,k    #numero do noh ij
        col4 = linha - Nz  # corresponde ao noh i,k-1  #numero do noh de cima

                
        c1 = +(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
        c4 = -(1.0/delta)
        
        RB = c1
        RBI = c4
        
        row[cont]  = linha
        col[cont]  = col1
        data[cont] = RB
        
        cont = cont + 1
            
#        row.append(linha)
#        col.append(col1)
#        data.append([np.complex128(RB)])
        
        
        row[cont]  = linha
        col[cont]  = col4
        data[cont] = RBI
        
        cont = cont + 1
            
#        row.append(linha)
#        col.append(col4)
#        data.append([np.complex128(RBI)])

    
    
    for i in range(1,Nx-1):
        
        # topo
        j = 0
        
        linha = (i)*Nz + j   #numero do noh ij
        
        col1 = linha	  # corresponde ao noh i,k    #numero do noh ij
        col3 = linha + 1  # corresponde ao noh i+1,k  #numero do noh da direita
        
        c1 = -(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
        c3 = +(1.0/delta)
        
        FS = c1
        FSI = c3

#                FS = 1.0
#                FSI = 1.0
        
#        row.append(linha)
#        col.append(col1)
#        data.append([np.complex128(FS)])
        
        
        row[cont]  = linha
        col[cont]  = col1
        data[cont] = FS
        
        cont = cont + 1
            
            
#        row.append(linha)
#        col.append(col3)
#        data.append([np.complex128(FSI)])
        
        row[cont]  = linha
        col[cont]  = col3
        data[cont] = FSI
        
        cont = cont + 1
            

        # fundo
        j = Nz - 1
        
        linha = (i)*Nz + j   #numero do noh ij
        
        col1 = linha	  # corresponde ao noh i,k    #numero do noh ij
        col2 = linha - 1  # corresponde ao noh i-1,k  #numero do noh da esquerda
                
        c1 = +(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
        c2 = -(1.0/delta)
        
        BB = c1
        BBI = c2


#        row.append(linha)
#        col.append(col1)
#        data.append([np.complex128(BB)])
        
        row[cont]  = linha
        col[cont]  = col1
        data[cont] = BB
        
        cont = cont + 1
            
        
#        row.append(linha)
#        col.append(col2)
#        data.append([np.complex128(BBI)])
        
        row[cont]  = linha
        col[cont]  = col2
        data[cont] = BBI
        
        cont = cont + 1
            
    
#    
    
    row = np.array(row)
#    print(row)
    
    col = np.array(col)
#    print(col)
    
    data = np.asarray(data)
    data = data.ravel()
#    print(data)
    
    
    Csp = csr_matrix((data, (row, col)))
#    Csp_ver = csr_matrix((data, (row, col)), shape=(Nx*Nz,Nx*Nz)).toarray()
#    print(Csp_ver)
    
    
    return Csp


