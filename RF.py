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

def Cal_N(isXrXl, isF, FS, lamda, fx, c, u, alpha, beta, phi, W, omega, kW, N):
    XlXr = [0] * Surface.shape
    tam1 = 0
    tam2 = 0
    tam3 = 0 
    tam4 = 0
    for i in range(Surface.shape):
        if isXrXl:
            tam1 = (c * beta[i] - u * beta[i] * math.tan(phi / 180 * math.pi)) * math.cos(alpha[i] / 180 * math.pi) / FS
            tam2 = N[i] * (math.tan(phi / 180 * math.pi) * math.cos(alpha[i] / 180 * math.pi) / FS - math.sin(alpha[i]) / 180 * math.pi)
            #       old N[i]s need to be add to list+params
            XlXr[i] = fx[i] * lamda * ( tam1 + tam2 - kW + D * math.cos(omega / 180 * math.pi))
        else:
            XlXr[i] = 0



def Cal_FS(isF):
    pass

def Calculating_FoS():
    pass 





#===========================================================================


 
