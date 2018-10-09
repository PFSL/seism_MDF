# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 19:54:25 2017

@author: Paulo
"""

import numpy as np


def solicita_vetor(Nx,Nz,delta,deep,num_rec,solicita,grade):

    """ ENTRADA DE DADOS AUTOMATICA DA RELACAO FONTES E RECEPTORES ###
    # 
    # ROTINA MONTAGEM DOS VETORES DE SOLICITACAO COM PULSO DE RICKER
    # COM SISTEMA DE AQUISIÇAO DO TIPO SPLIT-SPREAD COM NUMERO DE
    # RECEPTORES DEFINIDO NA ENTRADA DE DADOS
    #
    # NX        - NUMERO DE PONTOS EM X
    # NZ        - NUMERO DE PONTOS EM Z
    # DEEP      - PROFUNDIDADE DOS SENSORES
    # NUM_REC   - NUMERO DE RECEPTORES
    # SOLICITA  - MATRIZ COM FREQUENCIAS E CORRESPONDENTES AMPLITUDES
    # GRADE     - GRID DO DOMINIO A SER AVALIADO
    #
    #*****************************************************************"""

        
    
    
    [num_freq, lixo] = np.shape(solicita)
    print('\n\nnum_freq ',num_freq)
    
    print('\nsolicita: ', solicita)
    
    xi_off = (1 + num_rec/2)*delta
    xf_off = (Nx - 2 - num_rec/2)*delta
                    
    print('coord do primeiro no fonte do offset ', xi_off)
    print('coord do ultimo no fonte do offset ', xf_off)
    
    # varre grade com informacoes das fontes
    for i in range(0,num_freq):
        
        dist = np.zeros((2,Nx*Nz))
        
        
        for j in range(0,Nx*Nz):
            
            dist[0,j] = grade[j,0]
            
            dist[1,j] = np.sqrt((grade[j,1] - xi_off)**2 + (grade[j,2] - deep)**2)
            
#        print(dist)
            
        
        # ordena para tomar o no mais proximos
        dist = np.transpose(dist)
        
        node_sorted = dist[dist[:,1].argsort(),]

        node_si = int(node_sorted[0,0])

        
        
        
        dist = np.zeros((2,Nx*Nz))
        
        for j in range(0,Nx*Nz):
            
            dist[0,j] = grade[j,0]
            
            dist[1,j] = np.sqrt((grade[j,1] - xf_off)**2 + (grade[j,2] - deep)**2)
            
        
        # ordena para tomar o no mais proximos
        dist = np.transpose(dist)
        
        node_sorted = dist[dist[:,1].argsort(),]

        node_sf = int(node_sorted[0,0])
        
    
        
    
    

    print('\nno inicial com fonte atribuida')
    print(node_si)
    
    print('\nno final com fonte atribuida')
    print(node_sf)
    
    print('frequencia Hz')
    print(solicita[i][0])
    print('frequencia rad')
    print(solicita[i][0]*2*np.pi)
    
    print('amplitude')
    
    print(solicita[i][1])
    
        
    # numero de fontes disparadas
    num_src = node_sf - node_si + 1
    print('numero de fontes disparadas', num_src)
    
    # alocar memoria para vetor
    s = np.zeros((Nx*Nz,num_src*num_freq),dtype=np.complex128)

    omega_vec = []
    
    index_fones = np.zeros((num_rec,num_src*num_freq),dtype=np.int)
    
    cont = 0
    
    # varre grade com informacoes das fontes
    for i in range(0,num_freq):
        
        omega = solicita[i][0]*(2*np.pi)
        AA = solicita[i][1]
        
#        t0 = (1./omega)/np.pi
#        TAU = omega*t0/2.
        
        for src in range(0,num_src):
            
            fones_min = node_si + src - num_rec/2
            fones_max = node_si + src + num_rec/2
            print('minimo e maximo das fontes - ', fones_min, fones_max)
            
            vante = np.array(np.linspace(fones_min,node_si + src - 1,node_si + src - fones_min),dtype=np.int)
            print('vante - ',vante)
            
            re = np.array(np.linspace(node_si + src + 1,fones_max, fones_max - (node_si + src)),dtype=np.int)
            print('re - ',re)
            
            
            print(np.reshape(np.concatenate((vante,re),axis=0),(num_rec,1)))
            
            
            index_fones[:,cont] = np.array(np.reshape(np.concatenate((vante,re),axis=0),(num_rec,1))[:,0],dtype=np.int)
            
            print('index_fones completo',index_fones[:,cont])
            
            omega_vec.append(omega)
            
#            # ricker wavelet extraido do paper do edmundo "Efficient numerical models 
#            # for the prediction of acoustic wave propagation in the vicinity of a 
#            # wedge coastal region"
#            t0 = (1./omega)/np.pi
#            TAU = omega*t0/2.
#            
#            # calculo do valor do pulso de ricker
#            s(node_s,1) = AA*(2*sqrt(pi)*t0*exp(-imag*omega*ts)*TAU^2*exp(-TAU^2));
#            
#            # caluclo do pulso de ricker baseado no procedimento descrito pela franciane
#            s[node_s,0] = AA*omega**2*np.exp(-omega**2/(4*np.pi**2*omega**2))/(2*np.sqrt(2)*np.pi**3*(omega**2)**(3/2))
#            
            # valor do pulso adotado para testes
            s[node_si + src, cont] = s[node_si + i, cont] + AA
            
            
            cont = cont + 1
            
  
        
        
    return s, omega_vec, num_src, index_fones

