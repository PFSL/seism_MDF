import numpy as np
# import scipy as sp

from input_data import *
from truss_dyn import *



def f_obj_lbfgsb(P, *args):
    nodes = args[0]
    elem = args[1]
    bc = args[2]
    ti = args[3]
    tf = args[4]
    dt = args[5]
    ngl = args[6]
    sdof = args[7]
    nnel = args[8]
    nt = args[9]
    veli = args[10]
    desi = args[11]
    acc = args[12]
    vel = args[13]
    des = args[14]

    print("verifica tamanho do vetor forca")
    print(P)
    print(np.shape(P))

    [l, c] = np.shape(des)

    #    print("tamanho do argumento des")
    #    print(l,c)

    P = np.reshape(P, (l, 1))

    print("novo P")
    print(P)
    #



    # [des_lbfgs, acc_lbfgs, vel_lbfgs, stress_lbfgs, ngl, nnel, nt, nnos, nelem] = truss_dyn(nodes, elem, bc, ti, tf, dt, newF);
    des_opt, acc_opt, vel_opt, stress_opt, ngl, nnel, nt, nnos, nelem = truss_dyn(nodes, elem, bc, ti, tf, dt, P, ngl,
                                                                                  sdof, nnel, nt, veli, desi)

    #    print("des")
    #    print(des)
    #
    #    print("des_opt")
    #    print(des_opt)

    ## função residuo ---> mais importante
    error = (acc - acc_opt) + (vel - vel_opt) + (des - des_opt)
    # fob = (des - des_lbfgs)*10000
    # fob = vel - vel_lbfgs
    # fob = (acc - acc_lbfgs)
    # fob = (acc - acc_lbfgs) + 100*(vel - vel_lbfgs) + 10000*(des - des_lbfgs)
    print("error")
    print(error)
    print("")
    print(sum(error))
    print("")
    print(sum(sum(error)) ** 2)

    #    return sum(error)**2
    #    return np.linalg.norm(error)
    return sum(sum(error)) ** 2