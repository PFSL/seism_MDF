# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 12:53:21 2017

@author: Paulo
"""

import numpy as np
import os


# arquivo de saida para cada tiro dado no modelo
# nome do arquivo segue o padrao
# sol_MDF_acoustic + num_tiro + .dat
#
def wr_report_v2(b,K,Nx,Nz,grade,solicita,p,delta,Npml,Cpml):

    tam = np.size(p)

    # toma o diretorio atual do programa para retornar ao fim da rotina
    raiz = os.getcwd()
    saida = 'output'
    os.chdir(saida)
    
    
    arq_saida = 'sol_MDF_acoustic' + '.dat'

    leitura = open(arq_saida, 'w')
    
    
    leitura.write("\t*** F  I  N  I  T  E    D  I  F  F  E  R  E  N  C  E    M  E  T  H O  D ***\n\n")
    leitura.write("\t              *** 2  D    A  C  O  U  S  T  I  C  S ***\n\n")

    leitura.write("\n               P A R A M E T E R S   O F   D O M A I N\n")

    leitura.write("\n     XMIN                     XMAX                  1     4     7    x   ")
    leitura.write("\n     ,---------------------------, ZMIN             -*-----*-----*--->   ")
    leitura.write("\n     |                           |                  2|    5|    8|       ")
    leitura.write("\n  NZ |        D O M A I N        |                   *-----*-----*       ")
    leitura.write("\n     |                           |                  3|    6|    9|       ")
    leitura.write("\n     '---------------------------' ZMAX              *-----*-----*       ")
    leitura.write("\n                   NX                                |                   ")
    leitura.write("\n                                                   z\\/              \n\n")

    leitura.write(" XMIN  =  %9.3f \t\t ZMIN  =  %9.3f\n" % (np.min(grade[:,1]), np.min(grade[:,2])))
    leitura.write(" XMAX  =  %9.3f \t\t ZMAX  =  %9.3f\n" % (np.max(grade[:,1]), np.max(grade[:,2])))

    leitura.write("  NX   =  %9.0f \t\t  NZ  =  %9.0f\n\n" % (Nx, Nz))

    leitura.write(" NUMBER OF NODES  =  %9.0f\n" % (Nx*Nz))
    leitura.write(" DELTA            =  %9.3f\n" % (delta))

    # leitura.write("\nC O N D I C O E S   D E   C O N T O R N O   P R E S C R I T A S\n\n")
    leitura.write("\n\n     P R E S C R I B E D   S O U R C E S\n\n")

    leitura.write("\tNpml  =  %12.0f\n" % (Npml))
    leitura.write("\tCpml  =  %12.4f\n\n" % (Cpml))

#    [num_nos, num_dif_sources] = np.size(p)

#    leitura.write("\tNUMBER OF DIFFERENT SOURCES  =  %9.0f\n\n", num_dif_sources)

    leitura.write("\n   COORD X    COORD Z    FREQ (Hz)     AMP\n\n")

#    [num_f, lixo] = np.size(solicita)

#    for i in range(0,num_f):

    cx = solicita[0,0]
    cz = solicita[0,1]
    om = solicita[0,2]
    amp = solicita[0,3]

    leitura.write(" %9.3f  %9.3f  %9.3f  %9.3f\n" % (cx, cz, om, amp))



    leitura.write("\n                  M E D I A    P R O P E R T I E S\n")
    leitura.write("\n       ID\t  COORD X\tCOORD Z\t\tRHO\t  VELOCITY\n\n")

    count = 0
    for i in range(0,Nx):
        for k in range(0,Nz):
            
            dens = 1./b[i,k]
            velo = np.sqrt(K[i,k]*b[i,k])

            leitura.write("  %9.0f  %12.3f  %12.3f  %12.3f  %12.3f\n" % (count, grade[count,1], grade[count,2], dens, velo))
            
            count = count + 1

    leitura.write("\n\n P R E S S U R E   F O R   E A C H   F R E Q U E N C Y\n")



#    for i in range(0,num_dif_sources):
#
#        leitura.write("\n FREQUENCY  =  %12.3f\n", solicita[2]/(2*np.pi))

    leitura.write(" N O D E  \t P (real)\t P (imag)\n")

    for j in range(0,tam):

        pr = np.real(p[j])
        pi = np.imag(p[j])

        leitura.write("%10.0f\t%12.4f\t%12.4f\n" % (j, pr, pi))

        


#    leitura.write("\n\n P R E S S U R E   F O R   F R E Q U E N C Y   C O M P O S I T I O N\n\n")
#
#    leitura.write(" N O D E  \t\t\t P (real)\t\t\t P (imag)\n")
#
#    for i in range(0,num_nos):
#        
#        p_res_r = 0
#        p_res_i = 0
#        
#        for j in range(0,num_dif_sources):
#
#            p_res_r = p_res_r + real(p[i,j])
#            p_res_i = p_res_i + imag(p[i,j])
#
#        
#
#        leitura.write("%10.0f\t%12.4f\t%12.4f\n", i, p_res_r, p_res_i)




      
    
    leitura.close()
    os.chdir(raiz)
    
    return



































def escreve_res(Nx,Nz,grade,p,rho,v,omega,solicita,Npml,Cpml,arquivo):
    
    
    # create meshgrid
    [tam, lixo] = np.shape(p)
    
#    print(tam)
    
    xx = np.zeros((tam,1))
    zz = np.zeros((tam,1))
#    cont = 0;
    
    xx = []
    zz = []
    
    for line in grade:

        xx.append(float(line[0]))
        zz.append(float(line[1]))

    xx = np.array(xx)
    zz = np.array(zz)
    
#    for i in range(0,Nx):
#        for j in range(0,Nz):
#        
#            xx[cont,0] = X(0,i)
#            zz[cont,0] = Z(j,0)
#            cont = cont + 1
            
        
    # toma o diretorio atual do programa para retornar ao fim da rotina
    raiz = os.getcwd()
    saida = 'output'
    os.chdir(saida)

#    # nome do arquivo
#    arq_vtk = 'sol_MDF_acoustic_' + str(Nx) + '_' + str(Nz) + '.vtk'
#    arquivo = 'sol_MDF_acoustic_' + str(Nx) + '_' + str(Nz) + '.dat'
    
    print(arquivo)
    
    leitura = open(arquivo, 'w')
    
    leitura.write('\t\t*** M  E  T  O  D  O    D  A  S    D  I  F  E  R  E  N  C  A  S    F  I  N  I  T  A  S ***\n\n')
    leitura.write('\t\t                           *** A  C  U  S  T  I  C  A    2  D ***\n\n')
    
    leitura.write('\nP R O P R I E D A D E S   M A T E R I A I S\n')
    leitura.write('\nD E N S I D A D E S\n')
    aux = 1
    for prop in rho:
        leitura.write('R H O %i  =  %10.3f\n' % (aux, prop))
        aux = aux + 1
        
    aux = 1
    leitura.write('\nV E L O C I D A D E S\n')
    for prop in v:
        leitura.write('V E L O C I D A D E %1.0f  =  %10.3f\n' % (int(aux),prop))
        aux = aux + 1
        
    leitura.write('\nC O N D I C O E S   D E   C O N T O R N O   P R E S C R I T A S\n\n')
    
    [num_f, lixo] = np.shape(solicita)
    cont = 0
    for i in range(0,num_f):
        cont = cont + 1
        leitura.write('\tF O N T E     =  %10.1i\n' % cont)
        leitura.write('\tC O O R D   X =  %10.3f\n' % solicita[i,0])
        leitura.write('\tC O O R D   Z =  %10.3f\n' % solicita[i,1])
        leitura.write('\tO M E G A     =  %10.3f\n' % solicita[i,2])
        leitura.write('\tt s           =  %10.3f\n' % solicita[i,3])
        leitura.write('\tA M P         =  %10.3f\n' % solicita[i,4])
        leitura.write('\tN p m l       =  %10.3f\n' % Npml)
        leitura.write('\tC p m l       =  %10.3f\n\n' % Cpml)
        
    
    
    leitura.write('\n\t\tP R E S S O E S   N A S   C O O R D E N A D A S\n\n')
    leitura.write('  \tX\t\t\t Z\t\t\t P (real)\t\tP (imag)\t\tP (abs)\n')
    
    for i in range(0,tam):
        
        cx = xx[i]
        cy = zz[i]
        
        pr = np.real(p[i,0])
        pi = np.imag(p[i,0])
        pa = np.abs(p[i,0])
    
        leitura.write('%12.3f\t%12.3f\t%12.3f\t%12.3f\t%12.3f\n' % (cx, cy, pr, pi, pa))
    
    
    
    leitura.close()
    os.chdir(raiz)
    
    return
    