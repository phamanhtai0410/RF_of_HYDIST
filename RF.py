#Declaration 
import numpy as np
import pandas as pd
import os 
import sys



#Class definition
class Surface:
    def __init__(self, h, alpha, beta, phi, u, c, x, aa, W, kW, D, A ):
        self.h = h
        self.alpha = alpha 
        self.beta = beta
        self.phi = phi
        self.u = u
        self.c = c
        self.x = x
        self.aa = aa
        self.W = W
        self.kW = kW
        self.D = D
        self.A = A

class Normal_Force:
    def __init__(self, N_f, N_m):
        self.N_f = N_f
        self.N_m = N_m



#----------------------------------------------------------------------#
                #################################
                #      Calculating Functions    #
                #################################

def Cal_N(isXrXl, isF, FS, lamda, ):
    pass