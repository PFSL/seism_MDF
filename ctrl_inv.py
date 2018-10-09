#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 13:13:29 2018

@author: paulolomeu
"""

import numpy as np
import scipy as sp
import time


import input_data_new_1 as inp

import solicita_vetor as vec

import kernel_matrix_full as kfl

import kernel_matrix_sparse_opt as ksp

import kernel_matrix_sparse_opt_jit as ksj

#from report import *
#
from plot_pressao_ac_2d import *
#
#from helmholtz_fct import *



t0 = time.time()

# chamada das entradas do problema
[Nx, Nz, delta, grade, v, rho, Npml, Cpml, solicita, deep, num_rec, num_freq] = inp.input_data_new_1()

t1 = time.time()

# montagem do vetor de solicitacao
s, omega_vec, num_src, index_fones = vec.solicita_vetor(Nx,Nz,delta,deep,num_rec,solicita,grade)

t2 = time.time()

print(s)

print('\n\nTempo de execucao\n')
print('input time                 - %10.5f s' %(t1 - t0))
print('montar vetor solicitaca    - %10.5f s' %(t2 - t1))


omega = omega_vec[0]


[aux1, aux2] = np.shape(s)


p = np.zeros((aux1,aux2),dtype=np.complex128)

pspt = np.zeros((aux1,aux2),dtype=np.complex128)

psptjit = np.zeros((aux1,aux2),dtype=np.complex128)


fones = np.zeros((num_rec,num_freq*num_src),dtype=np.complex128)


t3 = time.time()

for i in range(0,num_freq*num_src):
    
    t4 = time.time()
#    print(omega_vec[i])
#    
#    
#    
#    print(max(s[:,i]))
#    
#    print(index_fones)
    
    
#    C = kfl.kernel_matrix_full(Nx,Nz,Npml,Cpml,omega,delta,rho,v)
    
    
    t5 = time.time()
    
    Cspt = ksp.kernel_matrix_sparse_opt(Nx,Nz,Npml,Cpml,omega,delta,rho,v)
    
#    Csptjit = ksj.kernel_matrix_sparse_opt_jit(Nx,Nz,Npml,Cpml,omega,delta,rho,v)
    
    t6 = time.time()
    
    # solver do sistema
#    p[:,i] = np.linalg.solve(C,s[:,i])
    
    
    t7 = time.time()
    
    pspt[:,i] = sp.sparse.linalg.spsolve(Cspt,s[:,i])
#    psptjit[:,i] = sp.sparse.linalg.spsolve(Csptjit,s[:,i])
    
    t8 = time.time()
    
    fones[:,i] = pspt[index_fones[:,i],i]
    
    
    print('montar matriz full      - %10.5f s' %(t5 - t4))
    
    print('montar matriz esparsa   - %10.5f s' %(t6 - t5))
    
    print('solver full system      - %10.5f s' %(t7 - t6))
    
    print('solver sparse system    - %10.5f s' %(t8 - t7))
    
    print('\n')
    ##print('solver system with gmres  - %10.5f s' %(t6 - t5))
    ##
    ##print('solver system with lgmres - %10.5f s' %(t7 - t6))
    ##
    ##print('solver system with SuperLU - %10.5f s' %(t8 - t7))
    
    
    
    
    
## nome do arquivo
#arq_vtk = 'sol_MDF_acoustic_' + str(Nx) + '.vtk'
##arq_res = 'sol_MDF_acoustic_' + str(Nx) + '_' + str(Nz) + '.dat'

#print(arq_vtk)

## plota pressoes
#plot_vtk_new(Nx, Nz, grade, pspt, arq_vtk, num_freq*num_src)


#print('\ndiferença dos resultados')
#print(np.linalg.norm(p-pspt))
#print(np.linalg.norm(psptjit-pspt))


#print(fones)


""""


COM A ETAPA ANTERIOR CONSEGUE MONTAR UM PROBLEMA GABARITO PARA SER RESOLVIDO
PELO SOLVER DO PROBLEMA INVERSO MAS O IDEAL AINDA É FAZER A MONTAGEM DA MATRIZ
PARA SOLUÇAO POR DIFERENÇAS FINITAS DE 5A ORDEM OU 13 PONTOS


FALTA ADAPTAR AGORA A ROTINA HELMHOLTZ PARA EXECUÇAO GERAR O RESULTADO
BASEADO NAS MEDIDAS DOS HIDROFONES, PARA ENTAO A ROTINA DE INVERSÃO SER
MODIFICADA COM OS PARAMETROS DO PROBLEMA EM QUESTAO


"""

#
#t1 = time.time()
#
## chamadas das rotinas de solucao
#t2 = time.time()
#
#Cspt = k_matrix_sp_t(Nx,Nz,Npml,Cpml,omega,b,K,delta,v)
#
#t3 = time.time()
#
## solver do sistema
#pspt = sp.sparse.linalg.spsolve(Cspt,s)
#
#t4 = time.time()
#
#
## nome do arquivo
#arq_vtk = 'sol_MDF_acoustic_' + str(Nx) + '_' + str(Nz) + '.vtk'
#arq_res = 'sol_MDF_acoustic_' + str(Nx) + '_' + str(Nz) + '.dat'
##
#
### outputs
### escreve arquivo de saida remover rho_stk e v_stk
### escreve_res(Nx,Nz,grade,p,rho_stk,v_stk,omega,solicita,Npml,Cpml,arq_res)
##wr_report_v2(b,K,Nx,Nz,grade,solicita,psp,delta,Npml,Cpml)
##
### # # verifica nos fonte e nos da grade adotados
### # verif_fontes(Nx,Nz,grade,solicita)
##
#
## plota pressoes
##plot_vtk(Nx, Nz, grade, p, arq_vtk)
#plot_vtk(Nx, Nz, grade, psp, arq_vtk)
#
##
##
##

##
#
#
#
#
#teste = helmholtz_2d(Nx, Nz, grade, b, K, Npml, Cpml, omega, delta, v, solicita)
#
#print(teste)
#
#
####
####
####print("")
####print("")
####print("")
####print("*********** INVERSAO COMECA AQUI ***********")
####
##### inversao para obtencao das forcas segundo um conjunto de medidas (acc, vel, des)
####P0 = np.ones((sdof, nt))
##### F0[2,:] = 350.
##### F0[3,:] = -1250.
##### F0[7,:] = -1000
##### F0[8,:] = 500
####
####
####args = (nodes, elem, bc, ti, tf, dt, ngl, sdof, nnel, nt, veli, desi, acc, vel, des)
####f_obj_lm(P0, *args)
####
####print("tamanho do P0")
####print(np.shape(P0))
####
##### l-BFGS-B
##### x, fob, d = optimize.fmin_l_bfgs_b(f_obj_lbfgsb, x0=F0, args=(nodes, elem, bc, ti, tf, dt, ngl, sdof, nnel, nt, veli, desi, acc, vel, des), approx_grad=True)
######print("vetor resultado")
######print(x)
##### print("valor da funcao objetivo")
##### print(fob)
##### print("flag resultado")
##### print(d)
####
####
##### sol = optimize.root(f_obj_lm, F0, args=(nodes, elem, bc, ti, tf, dt, ngl, sdof, nnel, nt, veli, desi, acc, vel, des), method='hybr', jac=False)
##### sol = optimize.root(f_obj_lm, F0, args=(nodes, elem, bc, ti, tf, dt, ngl, sdof, nnel, nt, veli, desi, acc, vel, des), method='lm', jac=False)
##### sol = optimize.root(f_obj_lm, F0, args=(nodes, elem, bc, ti, tf, dt, ngl, sdof, nnel, nt, veli, desi, acc, vel, des), method='broyden1', jac=False)
##### sol = optimize.root(f_obj_lm, F0, args=(nodes, elem, bc, ti, tf, dt, ngl, sdof, nnel, nt, veli, desi, acc, vel, des), method='broyden2', jac=False)
##### sol = optimize.root(f_obj_lm, F0, args=(nodes, elem, bc, ti, tf, dt, ngl, sdof, nnel, nt, veli, desi, acc, vel, des), method='anderson', jac=False)
####sol = optimize.root(f_obj_lm, F0, args=(nodes, elem, bc, ti, tf, dt, ngl, sdof, nnel, nt, veli, desi, acc, vel, des),
####                       method='krylov', jac=False)
####
####
##### import scipy
#####
##### scipy.optmi
####
####[l, c] = np.shape(des)
####
####xx = sol.x
####
####xx = np.reshape(sol.x, (l, c))
####print("vetor de forcas no tempo")
####print(xx)
####
####print("numero de iteraçoes")
####print(sol.nit)
####
####print("sucesso ou nao")
####print(sol.success)
####
####print("mensagem")
####print(sol.message)
####
####
####wr_report_inv(prob_name, xx)
####
####print("")
####print("*********** INVERSAO TERMINA AQUI ***********")
####
####
####
#
#
#




