#Declaration 
import numpy as np
import pandas as pd
import os 
import sys
import math


#Class definition

class Normal_Force:
    def __init__(self, N_f, N_m):
        self.N_f = N_f
        self.N_m = N_m


#======================================================================#
#                      PRE-PROCESSING 
########################################################################



########################################################################
#               1.        Class definition
########################################################################

class soil_props:
    def __init__(self, c, phi):
        self.c = c                  # Luc lien ket cua dat c 
        self.phi = phi              # Goc ma sat trong phi

class calculating_props:
    def __init__(self, u, kW, A, D):
        self.u = u
        self.kW = kW
        self.A = A
        self.D = D

class slide_props:
    def __init__(self, alpha, beta, a, x, W):
        self.alpha = alpha
        self.beta = beta
        self.a = a
        self.x = x
        self.W = W


        
#######################################################################
#                 2.        Pre-Processor Functions
#######################################################################

#----------------------------------------------------------------------

def parse_input(filename):
    #       Parsing matrix txt file to numpy array
    mat = np.loadtxt(filename)
    mat = np.transpose(mat)
    return mat

#----------------------------------------------------------------------




#----------------------------------------------------------------------
def index_first(cross_section_2D):
#       Defining index range of cross-section
    i = 0
    while cross_section_2D[0, i] == 0:
        i+=1
    first = i + 1 
    return first

#----------------------------------------------------------------------

def depth_converter(cross_section_3D):
#           Convert x y z cross-section type to x z cross section for RF pre- processor
    for i in range(cross_section_3D.shape[1]):
        cross_section_2D.append([math.sqrt((cross_section_3D[0, i] - cross_section_3D[0, 0]) ** 2 + (c[1, i] - cross_section_3D[1, 0]) ** 2) cross_section_3D[2, i]])
    return cross_section_2D
#----------------------------------------------------------------------

def center_defining(cross_section_2D, n):
#           Defining matrix of center points based on cross-section input
    center_list = [0 , 0] * n
    #first = index_first(cross_section_2D)
    

    ###



    ###

    if n ==1 :
        center_list[0, 0] = 
        center_list[1, 0] = 

    return center_list

#-----------------------------------------------------------------------

def entry_exit_point(cross_section_2D, center):
#           Find out the intersection points between cross-section and every circle
    R = Radius()
    i = 0


    return entry_point,exit_point

#------------------------------------------------------------------------

def radius_lines_defining(cross_section_2D, n, entry_RL, padding_RL):
#           Defining muti-radius_lines 
    radius_lines = [0] * n
    for i in range(n):
        radius_lines[i] = entry_RL + (i -1) * padding_RL
    return radius_lines

#-----------------------------------------------------------------------
#           Defining radius of every centers corresponding with every radius-lines
def Radius(center_y, RL):
    return R = abs(center_y - RL)

#-----------------------------------------------------------------------

#----------------------------------------------------------------------#
                #################################
                #      Calculating Functions    #
                #################################

def Cal_N(isXrXl, FS, lamda, fx, c, u, alpha, beta, phi, W, omega, kW):
    XlXr = [0] * Surface.shape
    tam1 = 0
    tam2 = 0
    tam3 = 0 
    tam4 = 0
    N = [0] * Surface.shape
    for i in range(Surface.shape):
        if isXrXl:
            tam1 = (c * beta[i] - u * beta[i] * math.tan(phi / 180 * math.pi)) * math.cos(alpha[i] / 180 * math.pi) / FS
            tam2 = N[i] * (math.tan(phi / 180 * math.pi) * math.cos(alpha[i] / 180 * math.pi) / FS - math.sin(alpha[i]) / 180 * math.pi)
            #       old N[i]s need to be add to list+params
            XlXr[i] = fx[i] * lamda * ( tam1 + tam2 - kW + D * math.cos(omega / 180 * math.pi))
        else:
            XlXr[i] = 0

        tam3 = (c * beta[i] * math.sin(alpha[i] / 180 * math.pi) + u * beta[i] * math.tan(phi / 180 * math.pi)) / FS
        tam4 = math.cos(alpha[i] / 180 * math.pi) + math.sin(alpha[i] / 180 * math.pi) + math.tan(phi / 180 * math.pi) / FS
        N[i] = (W[i] + XlXr[i] - tam3) / tam4
    return N
     
#--------------------------------------------------------------------------                    

def Cal_FS(isF, N, c, u, beta, alpha, phi, x, aa, A, R, W):
    tam1 = 0
    tam2 = 0
    tam3 = 0
    for i in range(Surface.shape):
        if isF:
            tam1 = tam1 + c * beta[i] * math.cos(alpha[i] / 180 * math.pi)
            tam2 = tam2 + (N[i] - u * beta[i] * math.tan(phi / 180 * math.pi) * math.cos(alpha[i] / 180 * math.pi))
            tam3 = tam3 + N[i] * math.sin(alpha[i] / 180 * math.pi) - A
        else
            tam1 = tam1 + tam1 + c * beta[i] * R
            tam2 = tam2 + (N[i] - u * beta[i]) * R * math.tan(phi / 180 * math.pi)
            tam3 = tam3 + W[i] * x[i] - A * aa[i]
    FS = (tam1 + tam2) / tam3
    return FS

#------------------------------------------------------------------------------

def Calculating_FoS(Surface, lamda, fx, Tolerance):
    
    # '////////////////////////////////////////////////////
    # '///////////////////////////////////////////////////
    # ' MAIN CALCULATING BELOW
    # ' Building an plot with every surface in vertical of lamda changing series and A list of value F_f,F_m
    # 'After got that, we take intersection point of two lines F_f(lamda) and F_m(lamda).
    # 'At that piont, we got F_m = F_f = F.O.S of the current surface
    # '   With every lamda values, We use RAPID SOLVER - Multi Iterations to have Convergence Status
    # '   The result of RAPID SOLVER is list of F_m corresponding with lamda and also with F_f
    # '


    # 'MAIN 4 STAGES of FOS CALCULATING
    # 'Main Concept of RAPID SOLVER - Newton Raphson Algorithm
    # '-----------------------------------------------------------------
    # '   STAGE 1: 1 Iteration : Set X = E = 0 at every slice ; set lamda = 0 ;
    # '   STAGE 2: 4-6 Iterations :    X = E = 0 as XlXr = False in Calculating_N ; to convergence F_m and F_f
    # '   STAGE 3: Main RAPID SOLVER - According to Newton Raphson Algorithm
    # '       set initial lamda = 2/3 slope between crest and toe
    # '       Selected Tolerance
    # '       Iterations until we got Abs(F_f - F_m) < Tolerance
    # '   STAGE 4: We got a list of F_f(lamda) and F_m(lamda)
    # '       Aferthat, we plot F_f(lamda) and F_m(lamda) ; the intersection is value of F.O.S that we finding
    
    #       DECLARATION AT THE BEGINNING POINT
    W[i] = [0] * Surface.shape
    N = [W[i] * alpha[i] for i in range(Surface.shape)]
    FS_f = [0] * lamda.shape
    FS_m = [0] * lamda.shape


    
    

    for i in range(lamda.shape):
        #   STAGE  1 : No X and E ; lamda = 0
        FS_f_crr = Cal_FS(True, N, c, u, beta, alpha, phi, x, aa, A. R, W)
        FS_m_crr = Cal_FS(False,N, c, u, beta, alpha, phi, x, aa, A. R, W)
        for j in range(fx.shape):
            #   STAGE 2
            for k in range(6):
                N = Cal_N(False, FS, lamda, fx, c, u, alpha, beta, phi, W, omega, kW)
                FS_f_crr = Cal_FS(True, N, c, u, beta, alpha, phi, x, aa, A, R, W)

            for k in range(6):
                N = Cal_N(False, FS, lamda, fx, c, u, alpha, beta, phi, W, omega, kW)
                FS_m_crr = Cal_FS(False, N, c, u, beta, alpha, phi, x, aa, A, R, W)

            #   STAGE 3
            while (abs(FS_f_crr - FS_m_crr) > Tolerance):
                N = Cal_N(True, FS_f_crr, lamda[i], fx[j], c, u, alpha, beta, phi, W, omega, kW)
                FS_f_crr = Cal_FS(True, N, c, u, beta, alpha, phi, x, aa, A, R, W)
                N = Cal_N(True, FS_m_crr, lamda[i], fx[j], c, u, alpha, beta, phi, W, omega, kW)
                FS_m_crr = Cal_FS(False, N, c, u, beta, alpha, phi, x, aa, A, R, W)
        FS_f[i] = FS_f_crr
        FS_m[i] = FS_m_crr

    #   STAGE 4: 
    


#===========================================================================


 
