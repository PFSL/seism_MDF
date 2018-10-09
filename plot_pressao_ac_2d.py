
# bibliotecas utilizadas
import numpy as np
import os


from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
from matplotlib.path import Path
import matplotlib.patches as patches









#===============================================================================
def plot_vtk_new(Nx, Nz, grade, pspt, arquivo, num_tiros):
    
    
    import numpy as np
    import os
    
#    import scipy
#    from scipy import ndimage
    
    # toma o diretorio atual do programa para retornar ao fim da rotina
    raiz = os.getcwd()
    saida = 'output'
    os.chdir(saida)
    
    print(arquivo)
    
    leitura = open(arquivo, 'w')
    
    leitura.write('# vtk DataFile Version 2.0\n')
        
    leitura.write('Structured Grid Example\n')
    leitura.write('ASCII\n')
    
    
#    leitura.write('DATASET STRUCTURED_POINTS\n')
    leitura.write('DATASET STRUCTURED_GRID\n')
#    leitura.write('DATASET RECTILINEAR_GRID\n')
    
    leitura.write('DIMENSIONS %i %i %i\n' % (Nx, 1, Nz))
    
#    X0 = 0.0
#    Y0 = 0.0
#    Z0 = 0.0
#    
#    leitura.write('ORIGIN %i %i %i\n' % (X0, Y0, Z0))
    nptos = Nx*Nz
    
    leitura.write('POINTS %i double\n' % nptos)
    
    for line in grade:
        
        leitura.write('%f %f %f\n' % (line[1],0,line[2]))
    
    
#    leitura.write('SPACING %f %f %f\n' % (DX, DY, DZ))
            
        
        
    
  
    
    leitura.write('POINT_DATA %i\n' % nptos)
    
    for k in range(0,num_tiros):
        
        leitura.write('SCALARS Pressao_tiro_' + str(k+1) + ' float\n' )
        leitura.write('LOOKUP_TABLE default\n')
        
        p = pspt[:,k]
        print(p)
        
        
        for line in p:
                
            # quando solver for esparso usar essa leitura
    #        print(line)
            pr = np.real(line)
            pi = np.imag(line)
            pa = np.abs(line)
                
            
    #        leitura.write('%12.8f' % pr)
            leitura.write('%14.11f\n' % pr)
    #        leitura.write('%10.3f %10.3f %10.3f\n' % (Cx, Cy, Cz))
        


    
    
    
    
    
#    
#    
#    leitura.write('DATASET STRUCTURED_POINTS\n')
##    leitura.write('DATASET STRUCTURED_GRID\n')
##    leitura.write('DATASET RECTILINEAR_GRID\n')
#    
#    leitura.write('DIMENSIONS %i %i %i\n' % (Nx, 1, Nz))
#    
#    X0 = 0.0
#    Y0 = 0.0
#    Z0 = 0.0
#    
#    leitura.write('ORIGIN %i %i %i\n' % (X0, Y0, Z0))
#    
#    
#    DX = (np.max(grade[:,1]) - np.min(grade[:,1]))/Nx
#    DY = 0
#    DZ = (np.max(grade[:,2]) - np.min(grade[:,2]))/Nz
#    
#    leitura.write('SPACING %f %f %f\n' % (DX, DY, DZ))
#            
#        
#        
#    nptos = Nx*Nz
#  
#    
#    leitura.write('POINT_DATA %i\n' % nptos)
#    
#    for k in range(0,num_tiros):
#        
#        leitura.write('SCALARS Pressao_tiro_' + str(k+1) + ' float 1\n' )
#        leitura.write('LOOKUP_TABLE default\n')
#        
#        p = pspt[:,k]
#        print(p)
#        
#        
#        for line in p:
#                
#            # quando solver for esparso usar essa leitura
#    #        print(line)
#            pr = np.real(line)
#            pi = np.imag(line)
#            pa = np.abs(line)
#                
#            
#    #        leitura.write('%12.8f' % pr)
#            leitura.write('%14.11f\n' % pr)
#    #        leitura.write('%10.3f %10.3f %10.3f\n' % (Cx, Cy, Cz))
#        
#        
##        for i in range(0,Nz):
##    
##            for j in range(0,Nx):
##       
##                line = (i)*Nz + j   #numero do noh ij
#                
#    #            # ordem das colunas, da menor p maior: col4, col2, col1, col3, col5
#    #            col1 = linha	  # corresponde ao noh i,k    #numero do noh ij
#    #            col2 = linha - 1  # corresponde ao noh i-1,k  #numero do noh da esquerda
#    #            col3 = linha + 1  # corresponde ao noh i+1,k  #numero do noh da direita
#    #            col4 = linha - Nz  # corresponde ao noh i,k-1  #numero do noh de cima
#    #            col5 = linha + Nz  # corresponde ao noh i,k+1  #numero do noh de baixo
#            
#            
##    #        for line in p:
##                
##                # quando solver for esparso usar essa leitura
##        #        print(line)
##                pr = np.real(p[line])
##                pi = np.imag(p[line])
##                pa = np.abs(p[line])
##                
##                
##        #        leitura.write('%12.8f' % pr)
##                leitura.write('%14.11f\n' % pr)
##        #        leitura.write('%10.3f %10.3f %10.3f\n' % (Cx, Cy, Cz))
#            
#        

    
    leitura.close()
    os.chdir(raiz)
    
    return












#===============================================================================
def plot_vtk(Nx, Nz, grade, p, arquivo):
    
    
    import numpy as np
    import os
    
#    import scipy
#    from scipy import ndimage
    
    # toma o diretorio atual do programa para retornar ao fim da rotina
    raiz = os.getcwd()
    saida = 'output'
    os.chdir(saida)
    
    
    leitura = open(arquivo, 'w')
    
    leitura.write('# vtk DataFile Version 2.0\n')
        
    leitura.write('Structured Grid Example\n')
    leitura.write('ASCII\n')
    leitura.write('DATASET STRUCTURED_POINTS\n')
#    leitura.write('DATASET STRUCTURED_GRID\n')
#    leitura.write('DATASET RECTILINEAR_GRID\n')
    
    leitura.write('DIMENSIONS %i %i %i\n' % (Nz, Nx, 1))
    
    X0 = 0.0
    Y0 = 0.0
    Z0 = 0.0
    
    leitura.write('ORIGIN %i %i %i\n' % (X0, Y0, Z0))
    
    
    DX = (np.max(grade[:,1]) - np.min(grade[:,1]))/Nx
    DY = (np.max(grade[:,2]) - np.min(grade[:,2]))/Nz
    DZ = 0
    
    leitura.write('SPACING %f %f %f\n' % (DX, DY, DZ))
            
        
        
    nptos = Nx*Nz
#    leitura.write('POINTS %i float\n' % nptos)
#    
#    for line in grade:
#        Cx = line[1]
#        Cy = line[2]
#        Cz = 0.0
#        leitura.write('%10.3f %10.3f %10.3f\n' % (Cx, Cy, Cz))
    
#    leitura.write('X_COORDINATES %i float\n' % nptos)
#    
#    for line in grade:
#        Cx = line[1]
#        leitura.write('%10.3f ' % Cx)
#    
#    leitura.write('\nY_COORDINATES %i float\n' % nptos)
#    
#    for line in grade:
#        Cy = line[2]
#        leitura.write('%10.3f ' % Cy)
#    
#    leitura.write('\nZ_COORDINATES %i float\n' % nptos)
#    
#    for line in grade:
#        Cz = 0.0
#        leitura.write('%10.3f ' % Cz)
    
    
    
    leitura.write('POINT_DATA %i\n' % nptos)
    leitura.write('SCALARS Pressao float 1\n')
    leitura.write('LOOKUP_TABLE default\n')
    
    for line in p:
#        # quando o solver e padrao usar essa leitura
#        pr = np.real(line[0])
#        pi = np.imag(line[0])
#        pa = np.abs(line[0])
        
        # quando solver for esparso usar essa leitura
#        print(line)
        pr = np.real(line)
        pi = np.imag(line)
        pa = np.abs(line)
        
        
#        leitura.write('%12.8f' % pr)
        leitura.write('%14.11f\n' % pr)
#        leitura.write('%10.3f %10.3f %10.3f\n' % (Cx, Cy, Cz))
    
    
    
    leitura.close()
    os.chdir(raiz)
    
    return












#===============================================================================
def le_arquivo(arquivo):
    
    # toma o diretorio atual do programa para retornar ao fim da rotina
    raiz = os.getcwd()
    saida = 'output'
    os.chdir(saida)
    
    # subrotina de leitura de arquivo
    leitura = open(arquivo, 'r')
    lista = []

    for line in leitura:
        lista.append(line)

    leitura.close()
    os.chdir(raiz)
    
    return lista


#===============================================================================
# rotina de leitura de dados
def read_inp(arq_input):
    
    # lista que armazena o arquivo de entrada (string)
    list_arq_ent = le_arquivo(arq_input)

    # indice referencia de cada informacao no arquivo de entrada (usar linha seguinte)
    ind_coord = list_arq_ent.index('		P R E S S O E S   N A S   C O O R D E N A D A S\n')
    
    coord = []
    # condicoes de cotorno atribuida para cada elemento do contorno (nao e para cada no)
    for line in list_arq_ent[ind_coord + 3:]:
        linha = line.split()
#        coord.append(map(float,linha))
        coord.append(linha)
        
#    print coord
    
    return coord
    

#===============================================================================
def plot_stress(grade, p):
    
    
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    import numpy as np
    import os
    
    import scipy
    from scipy import ndimage
    
    # toma o diretorio atual do programa para retornar ao fim da rotina
#    raiz = os.getcwd()
#    saida = 'output'
#    os.chdir(saida)
    
    
    
#    [num_f, aux_cr] = np.shape(list(coord))
#    
#    print(num_f, aux_cr)
#    
#    x = []
#    y = []
#    
#    for line in coord:
#
#        x.append(float(line[0]))
#        y.append(float(line[1]))
    
    x = grade[:,1]
    y = grade[:,2]
    
    print(x)
    print(y)
    
    print(p)
    pr = np.real(p)
    pi = np.imag(p)
    pa = np.abs(p)
    
    print(pr)
    print(pi)
    print(pa)
    
#    for line in coord:
#        pr.append(float(line[2]))
#        pi.append(float(line[3]))
#        pa.append(float(line[4]))


#    x.append(x[0])
#    y.append(y[0])
#    pr.append(pr[0])
#    pi.append(pi[0])
#    pa.append(pa[0])

#
#    x = np.transpose(x)
#    y = np.transpose(y)
#    
#    tam = np.shape(x)[0]
#    print('tam')
#    print(tam)
#    
#    pr = np.transpose(pr)
#    pi = np.transpose(pi)
#    pa = np.transpose(pa)
#    
##    x = np.array(np.reshape(x, (tam, 1)))
##    y = np.array(np.reshape(y, (tam, 1)))
#    
#    xx = x.flatten()
#    yy = y.flatten()
#    
#    pr = np.reshape(pr, (tam, 1))
#    pi = np.reshape(pi, (tam, 1))
#    pa = np.reshape(pa, (tam, 1))
#    
#    print(xx)
#    print(type(xx))
#    print(np.shape(xx))
#    print(yy)
#    print(type(yy))
#    print(np.shape(yy))
##    print(y)
##    print(pr)
##    print(pi)
##    print(pa)
#    
#    
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    triang = tri.Triangulation(x, y)
#    
#    
##    # tripcolor plot.
##    plt.figure()
##    plt.gca().set_aspect('equal')
##    plt.tripcolor(triang, pa, shading='flat', cmap=plt.cm.rainbow)
###    ndimage.rotate(figure, 90)
##    plt.colorbar()
##    plt.title('tripcolor of Delaunay triangulation, flat shading')
#    
    sizx = 10
    sizy = 5
    
    # Illustrate Gouraud shading.
    plt.figure(figsize=(sizx,sizy))
    plt.gca().set_aspect('equal')
    plt.tripcolor(triang, pa, shading='gouraud', cmap=plt.cm.rainbow)
    plt.colorbar()
    plt.title('PRESSAO ABSOLUTA')
    plt.savefig('abs1.png')
#        
#        
#    # pcolor plot.
#    plt.figure(figsize=(sizx,sizy))
#    plt.gca().set_aspect('equal')
#    plt.tricontourf(triang, pa)
##    plt.triplot(triang, color='0.97')
#    plt.colorbar()
#    plt.tricontour(triang, pa, colors='k')
#    plt.title('PRESSAO ABSOLUTA')
#    plt.savefig('abs2.png')
#    
#    
#    # Illustrate Gouraud shading.
#    plt.figure(figsize=(sizx,sizy))
#    plt.gca().set_aspect('equal')
#    plt.tripcolor(triang, pr, shading='gouraud', cmap=plt.cm.rainbow)
#    plt.colorbar()
#    plt.title('PRESSAO REAL')
#    plt.savefig('real1.png')
#        
#        
#    # pcolor plot.
#    plt.figure(figsize=(sizx,sizy))
#    plt.gca().set_aspect('equal')
#    plt.tricontourf(triang, pr)
##    plt.triplot(triang, color='0.97')
#    plt.colorbar()
#    plt.tricontour(triang, pr, colors='k')
#    plt.title('PRESSAO REAL')
#    plt.savefig('real2.png')
#    
#    
#    # Illustrate Gouraud shading.
#    plt.figure(figsize=(sizx,sizy))
#    plt.gca().set_aspect('equal')
#    plt.tripcolor(triang, pi, shading='gouraud', cmap=plt.cm.rainbow)
#    plt.colorbar()
#    plt.title('PRESSAO IMAGINARIA')
#    plt.savefig('imag1.png')
#        
#    # pcolor plot.
#    plt.figure(figsize=(sizx,sizy))
#    plt.gca().set_aspect('equal')
#    plt.tricontourf(triang, pi)
##    plt.triplot(triang, color='0.97')
#    plt.colorbar()
#    plt.tricontour(triang, pi, colors='k')
#    plt.title('PRESSAO IMAGINARIA')
#    plt.savefig('imag2.png')
#    
    
    
#    os.chdir(raiz)
    
    return