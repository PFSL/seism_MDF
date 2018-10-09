#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 20:35:01 2018

@author: paulolomeu
"""

import numpy as np


def kernel_check(Nx, Nz, Npml, Cpml, omega, b, K, delta):

    
    C = np.zeros((Nx*Nz,Nx*Nz))
 
    
    for i in range(0,Nx):
    
        for j in range(0,Nz):
    
            linha = (i)*Nz + j   #numero do noh ij
    
            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
            col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
            col2 = (i)*Nz + (j-1)  # corresponde ao noh i-1,j #numero do noh da esquerda
            col3 = (i)*Nz + (j+1)  # corresponde ao noh i+1,j	#numero do noh da direita
            col4 = ((i-1))*Nz + j # corresponde ao noh i,j-1  #numero do noh de cima
            col5 = ((i+1))*Nz + j  # corresponde ao noh i,j+1  #numero do noh de baixo
            

                        
            if(i != 0 and i != Nx-1 and j != 0 and j != Nz-1 ):
                C[linha,col1] = 1
                C[linha,col2] = 2
                C[linha,col3] = 3
                C[linha,col4] = 4
                C[linha,col5] = 5
                
                
            elif(i == 0):
                C[linha,col1] = 1
    #            C[linha,col2] = 2
                C[linha,col3] = 3
    #            C[linha,col4] = 4
    #            C[linha,col5] = 5
                

                
            elif(i == Nx-1):
                C[linha,col1] = 1
                C[linha,col2] = 2
    #            C[linha,col3] = 3
    #            C[linha,col4] = 4
    #            C[linha,col5] = 5
                
                
    
            elif(j == 0):
                C[linha,col1] = 1
    #            C[linha,col2] = 2
    #            C[linha,col3] = 3
    #            C[linha,col4] = 4
                C[linha,col5] = 5
                

                
            elif(j == Nz-1):
                C[linha,col1] = 1
    #            C[linha,col2] = 2
    #            C[linha,col3] = 3
                C[linha,col4] = 4
    #            C[linha,col5] = 5
                
           
    
    return C