# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 18:11:35 2017

@author: Paulo
"""

import numpy as np






def k_matrix_sp_t_old(Nx, Nz, Npml, Cpml, omega, b, K, delta, v):


    imag = 1j


    csix = np.zeros((Nx,2))
    csix = csix.view(dtype=np.complex128)


    for i in range(0,Nx):

        if(Npml > 1  and  i < Npml ): # dentro da primeira camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (i-1)/(Npml-1))

        elif(Npml > 1  and  i > (Nx-Npml+1)): # dentro da segunda camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (Nx-i)/(Npml-1))

        else: # dentro do modelo
          gammax = 0.0


        csix[i,0] = 1.0  + imag * gammax/omega

    #print csix


    csiz = np.zeros((Nx,2))
    csiz = csiz.view(dtype=np.complex128)

    for j in range(0,Nz):

        if(Npml > 1  and  j > (Nz-Npml+1)): # dentro da segunda camada de absorcao horizontal
          gammaz = Cpml*np.cos(np.pi/2.0 * np.float(Nz-j)/np.float(Npml-1))

        else: # dentro do modelo
          gammaz = 0.0


        csiz[j,0] = 1.0  + imag * gammaz/omega



    from scipy.sparse import csr_matrix    
    row = []
    col = []
    data = []


    for i in range(0,Nx):

        for j in range(0,Nz):

            linha = (i)*Nz + j   #numero do noh ij

            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
            col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
            col2 = (i)*Nz + (j-1)  # corresponde ao noh i-1,j #numero do noh da esquerda
            col3 = (i)*Nz + (j+1)  # corresponde ao noh i+1,j	#numero do noh da direita
            col4 = ((i-1))*Nz + j # corresponde ao noh i,j-1  #numero do noh de cima
            col5 = ((i+1))*Nz + j  # corresponde ao noh i,j+1  #numero do noh de baixo


            # no interno - nao tem c.c. aplicada
            if(i != 0 and i != Nx-1 and j != 0 and j != Nz-1 ):

                c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
                c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
                c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
                c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
                c1 = omega**2/K[i,j] - c2 -c3 - c4 - c5
                

                row.append(linha)
                col.append(col1)
                data.append(np.complex128(c1))
                
                row.append(linha)
                col.append(col2)
                data.append(np.complex128(c2))
                
                row.append(linha)
                col.append(col3)
                data.append(np.complex128(c3))
                
                row.append(linha)
                col.append(col4)
                data.append(np.complex128(c4))
                
                row.append(linha)
                col.append(col5)
                data.append(np.complex128(c5))
                 
                
                

            # borda esquerda
            elif(i == 0):

                LB = -v[i,j]/delta + omega*1j
                LBI = v[i,j]/delta


                row.append(linha)
                col.append(col1)
                data.append([np.complex128(LB)])
                
                row.append(linha)
                col.append(col3)
                data.append([np.complex128(LBI)])


            # borda direita
            elif(i == Nx-1):

                RB = v[i,j]/delta + omega*1j
                RBI = -v[i,j]/delta


                row.append(linha)
                col.append(col1)
                data.append([np.complex128(RB)])
                
                row.append(linha)
                col.append(col2)
                data.append([np.complex128(RBI)])


            # topo
            elif(j == 0):
                
#                FS = 1.0
#                FSI = -1.0
#                
                FS = v[i,j]/delta + omega*1j
                FSI = v[i,j]/delta
                
                row.append(linha)
                col.append(col1)
                data.append([np.complex128(FS)])
                
                row.append(linha)
                col.append(col5)
                data.append([np.complex128(FSI)])
                



            # fundo
            elif(j == Nz-1):

                BB = v[i,j]/delta + omega*1j
                BBI = -v[i,j]/delta


                row.append(linha)
                col.append(col1)
                data.append([np.complex128(BB)])
                
                row.append(linha)
                col.append(col4)
                data.append([np.complex128(BBI)])
    
    
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










# modelo seguindo arquivos do julia
# modificacoes em csix e csiz
# modificacoes em cond de contorno
def k_matrix_sp(Nx, Nz, Npml, Cpml, omega, b, K, delta, v):
    k_matrix_sp(Nx, Nz, Npml, Cpml, omega, b, K, delta, v)


    imag = 1j


    csix = np.zeros((Nx,1), dtype=np.complex128)
#    csix = csix.view(dtype=np.complex128)


    for i in range(0,Nx):

        if(Npml > 1  and  i < Npml ): # dentro da primeira camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (i-1)/(Npml-1))

        elif(Npml > 1  and  i > (Nx-Npml)): # dentro da segunda camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (Nx-i)/(Npml-1))

        else: # dentro do modelo
          gammax = 0.0


        csix[i] = 1.0  + imag * gammax/omega



    csiz = np.zeros((Nx,1), dtype=np.complex128)
#    csiz = csiz.view(dtype=np.complex128)

    for k in range(0,Nz):
        
        
#        if(Npml > 1  and  j < Npml ): # dentro da primeira camada de absorcao vertical (topo)
#          gammax = Cpml*np.cos(np.pi/2.0 * (j-1)/(Npml-1))
#
#        elif(Npml > 1  and  j > (Nx-Npml+1)): # dentro da segunda camada de absorcao vertical (leito)
#          gammax = Cpml*np.cos(np.pi/2.0 * (Nx-j)/(Npml-1))
#          

        # so tem um else se so tiver a camada inferior
        if (Npml > 1  and  k > (Nz-Npml)): # dentro da segunda camada de absorcao horizontal (leito)
            gammaz = Cpml*np.cos(np.pi/2.0 * (Nz-k)/(Npml-1))

        elif (Npml > 1  and  k < Npml ): # dentro da primeira camada de absorcao vertical (topo)
            gammaz = Cpml*np.cos(np.pi/2.0 * (k-1)/(Npml-1))

        else: # dentro do modelo
            gammaz = 0.0
          


        csiz[k] = 1.0  + imag * gammaz/omega



    from scipy.sparse import csr_matrix    
    row = []
    col = []
    data = []


    for i in range(0,Nx):

        for j in range(0,Nz):

            linha = (i)*Nz + j   #numero do noh ij

#            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
#            col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
#            col2 = (i)*Nz + (j-1)  # corresponde ao noh i-1,j #numero do noh da esquerda
#            col3 = (i)*Nz + (j+1)  # corresponde ao noh i+1,j	#numero do noh da direita
#            col4 = ((i-1))*Nz + j # corresponde ao noh i,j-1  #numero do noh de cima
#            col5 = ((i+1))*Nz + j  # corresponde ao noh i,j+1  #numero do noh de baixo
            
            
            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
            col1 = linha	  # corresponde ao noh i,k    #numero do noh ij
            col2 = linha - 1  # corresponde ao noh i-1,k  #numero do noh da esquerda
            col3 = linha + 1  # corresponde ao noh i+1,k  #numero do noh da direita
            col4 = linha - Nz  # corresponde ao noh i,k-1  #numero do noh de cima
            col5 = linha + Nz  # corresponde ao noh i,k+1  #numero do noh de baixo



            # no interno - nao tem c.c. aplicada
            if(i != 0 and i != Nx-1 and j != 0 and j != Nz-1 ):
                # staggered grid five point stencil (hustedt et al., 2004 equacao 8)
                c2 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i-1,j] )/(csiz[j]+csiz[j-1])
                c3 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i+1,j] )/(csiz[j]+csiz[j+1])
                c4 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i,j-1] )/(csix[i]+csix[i-1])
                c5 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i,j+1] )/(csix[i]+csix[i+1])
                c1 = omega**2/K[i,j] - c2 -c3 - c4 - c5

#                c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
#                c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
#                c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
#                c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
#                c1 = omega**2/K[i,j] - c2 -c3 - c4 - c5
                

                row.append(linha)
                col.append(col1)
                data.append(np.complex128(c1))
                
                row.append(linha)
                col.append(col2)
                data.append(np.complex128(c2))
                
                row.append(linha)
                col.append(col3)
                data.append(np.complex128(c3))
                
                row.append(linha)
                col.append(col4)
                data.append(np.complex128(c4))
                
                row.append(linha)
                col.append(col5)
                data.append(np.complex128(c5))
                 
                
                

            # borda esquerda
            elif(i == 0):
                
                c1 = -(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
                c5 = +(1.0/delta)
                
                LB = c1
                LBI = c5


                row.append(linha)
                col.append(col1)
                data.append([np.complex128(LB)])
                
                row.append(linha)
                col.append(col5)
                data.append([np.complex128(LBI)])


            # borda direita
            elif(i == Nx-1):
                
                c1 = +(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
                c4 = -(1.0/delta)
                
                RB = c1
                RBI = c4

#                RB = v[i,j]/delta + omega*1j
#                RBI = -v[i,j]/delta

                row.append(linha)
                col.append(col1)
                data.append([np.complex128(RB)])
                
                row.append(linha)
                col.append(col4)
                data.append([np.complex128(RBI)])


            # topo
            elif(j == 0):
                
                c1 = -(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
                c3 = +(1.0/delta)
                
                FS = c1
                FSI = c3

#                FS = 1.0
#                FSI = 1.0
                
                row.append(linha)
                col.append(col1)
                data.append([np.complex128(FS)])
                
                row.append(linha)
                col.append(col3)
                data.append([np.complex128(FSI)])
                


            # fundo
            elif(j == Nz-1):
                
                c1 = +(1.0/delta + 1j * omega/np.sqrt(K[i,j]*b[i,j]))
                c2 = -(1.0/delta)
                
                BB = c1
                BBI = c2

#                BB = v[i,j]/delta + omega*1j
#                BBI = v[i,j]/delta

                row.append(linha)
                col.append(col1)
                data.append([np.complex128(BB)])
                
                row.append(linha)
                col.append(col2)
                data.append([np.complex128(BBI)])
    
    
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

























# acertar os valores de cond contorno modificadas para teste
def k_matrix_sp_old(Nx, Nz, Npml, Cpml, omega, b, K, delta, v):


    imag = 1j


    csix = np.zeros((Nx,1), dtype=np.complex128)
    csix = csix.view(dtype=np.complex128)


    for i in range(0,Nx):

        if(Npml > 1  and  i < Npml ): # dentro da primeira camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (i-1)/(Npml-1))

        elif(Npml > 1  and  i > (Nx-Npml+1)): # dentro da segunda camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (Nx-i)/(Npml-1))

        else: # dentro do modelo
          gammax = 0.0


        csix[i] = 1.0  + imag * gammax/omega

    #print csix


    csiz = np.zeros((Nx,1), dtype=np.complex128)
    csiz = csiz.view(dtype=np.complex128)

    for j in range(0,Nz):

        if(Npml > 1  and  j > (Nz-Npml+1)): # dentro da segunda camada de absorcao horizontal
          gammaz = Cpml*np.cos(np.pi/2.0 * np.float(Nz-j)/np.float(Npml-1))

        else: # dentro do modelo
          gammaz = 0.0


        csiz[j] = 1.0  + imag * gammaz/omega



    from scipy.sparse import csr_matrix    
    row = []
    col = []
    data = []


    for i in range(0,Nx):

        for j in range(0,Nz):

            linha = (i)*Nz + j   #numero do noh ij

            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
            col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
            col2 = (i)*Nz + (j-1)  # corresponde ao noh i-1,j #numero do noh da esquerda
            col3 = (i)*Nz + (j+1)  # corresponde ao noh i+1,j	#numero do noh da direita
            col4 = ((i-1))*Nz + j # corresponde ao noh i,j-1  #numero do noh de cima
            col5 = ((i+1))*Nz + j  # corresponde ao noh i,j+1  #numero do noh de baixo


            # no interno - nao tem c.c. aplicada
            if(i != 0 and i != Nx-1 and j != 0 and j != Nz-1 ):

                c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
                c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
                c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
                c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
                c1 = omega**2/K[i,j] - c2 -c3 - c4 - c5
                

                row.append(linha)
                col.append(col1)
                data.append(np.complex128(c1))
                
                row.append(linha)
                col.append(col2)
                data.append(np.complex128(c2))
                
                row.append(linha)
                col.append(col3)
                data.append(np.complex128(c3))
                
                row.append(linha)
                col.append(col4)
                data.append(np.complex128(c4))
                
                row.append(linha)
                col.append(col5)
                data.append(np.complex128(c5))
                 
                
                

            # borda esquerda
            elif(i == 0):
                
#                c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
                c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
#                c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
#                c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
                c1 = omega**2/K[i,j] - c3
                
                print(c1,c3)
                
#                LB = c1
                LB = 1
#                LBI = c3
                LBI = 1

#                LB = -v[i,j]/delta + omega*1j
#                LBI = v[i,j]/delta
                print(LB)


                row.append(linha)
                col.append(col1)
                data.append([np.complex128(LB)])
                
                row.append(linha)
                col.append(col3)
                data.append([np.complex128(LBI)])


            # borda direita
            elif(i == Nx-1):
                
                c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
#                c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
#                c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
#                c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
                c1 = omega**2/K[i,j] - c2
                
#                RB = c1
                RB = 1
#                RBI = c2
                RBI = 1

#                RB = v[i,j]/delta + omega*1j
#                RBI = -v[i,j]/delta

                row.append(linha)
                col.append(col1)
                data.append([np.complex128(RB)])
                
                row.append(linha)
                col.append(col2)
                data.append([np.complex128(RBI)])


            # topo
            elif(j == 0):
                
#                c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
#                c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
#                c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
                c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
                c1 = omega**2/K[i,j] - c5
                
#                FS = c1
                FS = 1
#                FSI = c5
                FSI = 1

#                FS = 1.0
#                FSI = 1.0
                
                row.append(linha)
                col.append(col1)
                data.append([np.complex128(FS)])
                
                row.append(linha)
                col.append(col5)
                data.append([np.complex128(FSI)])
                


            # fundo
            elif(j == Nz-1):
                
#                c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
#                c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
                c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
#                c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
                c1 = omega**2/K[i,j] - c4
                
#                BB = c1
                BB = 1
#                BBI = c4
                BBI = 10

#                BB = v[i,j]/delta + omega*1j
#                BBI = v[i,j]/delta

                row.append(linha)
                col.append(col1)
                data.append([np.complex128(BB)])
                
                row.append(linha)
                col.append(col4)
                data.append([np.complex128(BBI)])
    
    
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











def k_sp_div(Nx, Nz, Npml, Cpml, omega, b, K, delta, v):


    imag = 1j


    csix = np.zeros((Nx,2))
    csix = csix.view(dtype=np.complex128)


    for i in range(0,Nx):

        if(Npml > 1  and  i < Npml ): # dentro da primeira camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (i-1)/(Npml-1))

        elif(Npml > 1  and  i > (Nx-Npml+1)): # dentro da segunda camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (Nx-i)/(Npml-1))

        else: # dentro do modelo
          gammax = 0.0


        csix[i,0] = 1.0  + imag * gammax/omega

    #print csix


    csiz = np.zeros((Nx,2))
    csiz = csiz.view(dtype=np.complex128)

    for j in range(0,Nz):

        if(Npml > 1  and  j > (Nz-Npml+1)): # dentro da segunda camada de absorcao horizontal
          gammaz = Cpml*np.cos(np.pi/2.0 * np.float(Nz-j)/np.float(Npml-1))

        else: # dentro do modelo
          gammaz = 0.0


        csiz[j,0] = 1.0  + imag * gammaz/omega



    from scipy.sparse import csr_matrix    
    row = []
    col = []
    data = []


    for i in range(1,Nx-1):

        for j in range(1,Nz-1):

            linha = (i)*Nz + j   #numero do noh ij

            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
            col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
            col2 = (i)*Nz + (j-1)  # corresponde ao noh i-1,j #numero do noh da esquerda
            col3 = (i)*Nz + (j+1)  # corresponde ao noh i+1,j	#numero do noh da direita
            col4 = ((i-1))*Nz + j # corresponde ao noh i,j-1  #numero do noh de cima
            col5 = ((i+1))*Nz + j  # corresponde ao noh i,j+1  #numero do noh de baixo


#            # no interno - nao tem c.c. aplicada
#            if(i != 0 and i != Nx-1 and j != 0 and j != Nz-1 ):

            c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
            c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
            c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
            c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
            c1 = omega**2/K[i,j] - c2 -c3 - c4 - c5
            

            row.append(linha)
            col.append(col1)
            data.append(np.complex128(c1))
            
            row.append(linha)
            col.append(col2)
            data.append(np.complex128(c2))
            
            row.append(linha)
            col.append(col3)
            data.append(np.complex128(c3))
            
            row.append(linha)
            col.append(col4)
            data.append(np.complex128(c4))
            
            row.append(linha)
            col.append(col5)
            data.append(np.complex128(c5))
                 
                
    for j in range(0,Nz):            
        
        i = 0 # borda esquerda
            
        linha = (i)*Nz + j   #numero do noh ij

        # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
        col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
        col3 = (i)*Nz + (j+1)  # corresponde ao noh i+1,j	#numero do noh da direita

        LB = -v[i,j]/delta + omega*1j
        LBI = v[i,j]/delta

        row.append(linha)
        col.append(col1)
        data.append([np.complex128(LB)])
        
        row.append(linha)
        col.append(col3)
        data.append([np.complex128(LBI)])



    for j in range(0,Nz):            
        
        i = Nx-1 # borda direita
        
        linha = (i)*Nz + j   #numero do noh ij

        # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
        col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
        col2 = (i)*Nz + (j-1)  # corresponde ao noh i-1,j #numero do noh da esquerda
        
        RB = v[i,j] + omega*1j
        RBI = -v[i,j]/delta

        row.append(linha)
        col.append(col1)
        data.append([np.complex128(RB)])
        
        row.append(linha)
        col.append(col2)
        data.append([np.complex128(RBI)])



    for i in range(1,Nx-1):
        
        j == 0 # topo
        
        linha = (i)*Nz + j   #numero do noh ij

        # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
        col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
        col5 = ((i+1))*Nz + j  # corresponde ao noh i,j+1  #numero do noh de baixo
        

#        FS = 1.0
#        row.append(linha)
#        col.append(col1)
#        data.append([np.complex128(FS)])
#        
#        row.append(linha)
#        col.append(col5)
#        data.append([np.complex128(FS)])
#            
        
        FS = -v[i,j] + omega*1j
        FSI = v[i,j]/delta
        
        row.append(linha)
        col.append(col1)
        data.append([np.complex128(FS)])
        
        row.append(linha)
        col.append(col5)
        data.append([np.complex128(FSI)])
            


    for i in range(1,Nx-1):
        
        j == Nz-1 # fundo
        
        linha = (i)*Nz + j   #numero do noh ij

        # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
        col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
        col4 = ((i-1))*Nz + j # corresponde ao noh i,j-1  #numero do noh de cima
        

        BB = v[i,j] + omega*1j
        BBI = -v[i,j]/delta


        row.append(linha)
        col.append(col1)
        data.append([np.complex128(BB)])
        
        row.append(linha)
        col.append(col4)
        data.append([np.complex128(BBI)])
    
    

    row = np.array(row)
#    print(row)
    
#    print('col')
#    print(col)
    col = np.array(col)
#    print(col)
    
#    print('data')
#    print(data)
    data = np.asarray(data)
#    data = np.array(data, dtype=np.complex128)
    data = data.ravel()
#    print(data)
    
    
    Csp = csr_matrix((data, (row, col)))
#    Csp_ver = csr_matrix((data, (row, col)), shape=(Nx*Nz,Nx*Nz)).toarray()
#    print(Csp_ver)
    
    
    return Csp













#
#import scipy as sp

#from memory_profiler import profile
#
#@profile
def k_matrix_new_sp(Nx, Nz, Npml, Cpml, omega, b, K, delta, v):


    imag = 1j


    csix = np.zeros((Nx,2))
    csix = csix.view(dtype=np.complex128)


    for i in range(0,Nx):

        if(Npml > 1  and  i < Npml ): # dentro da primeira camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (i-1)/(Npml-1))

        elif(Npml > 1  and  i > (Nx-Npml+1)): # dentro da segunda camada de absorcao vertical
          gammax = Cpml*np.cos(np.pi/2.0 * (Nx-i)/(Npml-1))

        else: # dentro do modelo
          gammax = 0.0


        csix[i,0] = 1.0  + imag * gammax/omega

    #print csix


    csiz = np.zeros((Nx,2))
    csiz = csiz.view(dtype=np.complex128)

    for j in range(0,Nz):

        if(Npml > 1  and  j > (Nz-Npml+1)): # dentro da segunda camada de absorcao horizontal
          gammaz = Cpml*np.cos(np.pi/2.0 * np.float(Nz-j)/np.float(Npml-1))

        else: # dentro do modelo
          gammaz = 0.0


        csiz[j,0] = 1.0  + imag * gammaz/omega

    #print csiz



    
    #C = np.zeros((Nx*Nz,Nx*Nz))
    C = np.zeros((Nx*Nz,2*Nx*Nz))
    C = C.view(dtype=np.complex128)
#    print(C)


    from scipy.sparse import csr_matrix    
#    Csp = coo_matrix((Nx*Nz, Nx*Nz), dtype=np.complex128)
    row = []
    col = []
    data = []


    for i in range(0,Nx):

        for j in range(0,Nz):

            linha = (i)*Nz + j   #numero do noh ij

            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
            col1 = (i)*Nz + j	   # corresponde ao noh i,j  #numero do noh ij
            col2 = (i)*Nz + (j-1)  # corresponde ao noh i-1,j #numero do noh da esquerda
            col3 = (i)*Nz + (j+1)  # corresponde ao noh i+1,j	#numero do noh da direita
            col4 = ((i-1))*Nz + j # corresponde ao noh i,j-1  #numero do noh de cima
            col5 = ((i+1))*Nz + j  # corresponde ao noh i,j+1  #numero do noh de baixo


            # no interno - nao tem c.c. aplicada
            if(i != 0 and i != Nx-1 and j != 0 and j != Nz-1 ):

                c2 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i-1,j] )/(csix[i]+csix[i-1])
                c3 = (1.0/(csix[i] * delta**2))  * (b[i,j]+b[i+1,j] )/(csix[i]+csix[i+1])
                c4 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j-1] )/(csiz[j]+csiz[j-1])
                c5 = (1.0/(csiz[j] * delta**2))  * (b[i,j]+b[i,j+1] )/(csiz[j]+csiz[j+1])
                c1 = omega**2/K[i,j] - c2 -c3 - c4 - c5
                
#                c1 = np.complex128(c1)
#                c2 = np.complex128(c2)
#                c3 = np.complex128(c3)
#                c4 = np.complex128(c4)
#                c5 = np.complex128(c5)

                C[linha,col1:col1+1] = c1
#                print(c1)
                row.append(linha)
                col.append(col1)
                data.append(c1)
                
                C[linha,col2:col2+1] = c2
#                print(c2)
                row.append(linha)
                col.append(col2)
                data.append(c2)
                
                C[linha,col3:col3+1] = c3
#                print(c3)
                row.append(linha)
                col.append(col3)
                data.append(c3)
                
                C[linha,col4:col4+1] = c4
#                print(c4)
                row.append(linha)
                col.append(col4)
                data.append(c4)
                
                C[linha,col5:col5+1] = c5
#                print(c5)
                row.append(linha)
                col.append(col5)
                data.append(c5)
                 
                
                

            # borda esquerda
            elif(i == 0):

                LB = -v[i,j]/delta + omega*1j
                LBI = v[i,j]/delta

                C[linha,col1:col1+1] = np.complex128(LB)
#                print([np.complex128(LB)])
                row.append(linha)
                col.append(col1)
                data.append([np.complex128(LB)])
                
                C[linha,col3:col3+1] = np.complex128(LBI)
#                print([np.complex128(LBI)])
                row.append(linha)
                col.append(col3)
                data.append([np.complex128(LBI)])


            # borda direita
            elif(i == Nx-1):

                RB = v[i,j] + omega*1j
                RBI = -v[i,j]/delta

                C[linha,col1:col1+1] = np.complex128(RB)
#                print([np.complex128(RB)])
                row.append(linha)
                col.append(col1)
                data.append([np.complex128(RB)])
                
                C[linha,col2:col2+1] = np.complex128(RBI)
#                print([np.complex128(RBI)])
                row.append(linha)
                col.append(col2)
                data.append([np.complex128(RBI)])


            # topo
            elif(j == 0):
                       
                FS = 1.0

                C[linha,col1:col1+1] = np.complex128(FS)
#                print([np.complex128(FS)])
                row.append(linha)
                col.append(col1)
                data.append([np.complex128(FS)])
                
                C[linha,col5:col5+1] = np.complex128(FS)
#                print([np.complex128(FS)])
                row.append(linha)
                col.append(col5)
                data.append([np.complex128(FS)])


            # fundo
            elif(j == Nz-1):

                BB = v[i,j] + omega*1j
                BBI = v[i,j]/delta

                C[linha,col1:col1+1] = np.complex128(BB)
#                print([np.complex128(BB)])
                row.append(linha)
                col.append(col1)
                data.append([np.complex128(BB)])
                
                C[linha,col4:col4+1] = np.complex128(BBI)
#                print([np.complex128(BBI)])
                row.append(linha)
                col.append(col4)
                data.append([np.complex128(BBI)])
    
    
#    Nval = np.size(row)
#    
#    dataarr = np.zeros((Nval), dtype=np.complex128)
#    cont = 0
#    
#    for item in data:
#        
#        dataarr[cont] = item
#        
#        cont = cont + 1
        
    
    
#    print('rows')
#    print(row)
    row = np.array(row)
#    print(row)
    
#    print('col')
#    print(col)
    col = np.array(col)
#    print(col)
    
#    print('data')
#    print(data)
    data = np.asarray(data)
#    data = np.array(data, dtype=np.complex128)
    data = data.ravel()
#    print(data)
    
    
    Csp = csr_matrix((data, (row, col)))
    Csp_ver = csr_matrix((data, (row, col)), shape=(Nx*Nz,Nx*Nz)).toarray()
#    print(Csp_ver)
    
    
    return C, Csp, Csp_ver



