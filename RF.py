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

class setting:
    def __init__(self, fx, lamda, Tolerance):
        self.fx = fx
        self.lamda = lamda
        self.Tolerance = Tolerance

class soil_props:
    def __init__(self, c, phi, gamma):
        self.c = c                  # Luc lien ket cua dat c 
        self.phi = phi              # Goc ma sat trong phi
        self.gamma = gamma          # Trong luong rieng cua dat

class calculating_props:
    def __init__(self, u, kW, A, D, omega):
        self.u = u
        self.kW = kW
        self.A = A
        self.D = D
        self.omega = omega

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
#       Cross-section is a numpy array
    
    if cross_section_2D[0, 1] == 0:
        isLeft = True
    else:
        isLeft = False
    if isLeft == False:
        cross_section_2D = np.flip(cross_section_2D, axis = 0)
    i = 0
    while cross_section_2D[i, 1] == 0:
        i+=1
    
    if isLeft:
        return isLeft, i
    else:
        return isLeft, (cross_section_2D.shape[0] - i - 1)


#----------------------------------------------------------------------
# Receive a list of surface depth 
# Return a numpy array of depth converter

def depth_converter(cross_section_3D):
    cross_section_3D = np.array(cross_section_3D)
    cross_section_2D = list()
    for i in range(cross_section_3D.shape[0]):
        temp = [math.sqrt((cross_section_3D[i, 0] - cross_section_3D[0, 0]) ** 2 + (cross_section_3D[i, 1] - cross_section_3D[0, 1])**2), cross_section_3D[i, 2]]
        cross_section_2D.append(temp)
    return np.array(cross_section_2D)

#----------------------------------------------------------------------
#  Receive a numpy array of depth
#  Return an array of center
def center_defining(cross_section_2D, n, first, isLeft, dx):
#           Defining matrix of center points based on cross-section input
    center_arr = np.array([0 , 0] * (n**2)).reshape(n**2, 2)
    note = None
    if isLeft:
        note = (cross_section_2D.shape[0] - first) // 2
    else:
        note = first // 2
    
    for i in range(n):
        for j in range(n):
            center_arr[i * n + j, 0] = cross_section_2D[note, 0] - (n // 2 - j) * dx
            center_arr[i * n + j, 1] = -i * dx

    return center_arr

#-----------------------------------------------------------------------

#------------------------------------------------------------------------

def radius_lines_defining(cross_section_2D, n, entry_RL, padding_RL):
#           Defining multi-radius_lines 
    radius_lines = np.array([0] * n)
    for i in range(n):
        radius_lines[i] = entry_RL + (i -1) * padding_RL
    return radius_lines
    
#-----------------------------------------------------------------------
#           Defining radius of every centers corresponding with every radius-lines
def Radius(center_y, RL):
    return abs(center_y - RL)


#-----------------------------------------------------------------------
def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y



#-----------------------------------------------------------------------

#----------------------------------------------------------------------#
                #################################
                #      Calculating Functions    #
                #################################

def Cal_N(Surface, isXrXl, FS, lamda, fx, c, u, alpha, beta, phi, W, omega, kW):
    XlXr = [0] * Surface.shape[0]
    tam1 = 0
    tam2 = 0
    tam3 = 0 
    tam4 = 0
    N = [0] * Surface.shape[0]
    for i in range(Surface.shape[0] - 1):
        if isXrXl:
            tam1 = (c * beta[i] - u * beta[i] * math.tan(phi)) * math.cos(alpha[i]) / FS
            tam2 = N[i] * (math.tan(phi) * math.cos(alpha[i] ) / FS - math.sin(alpha[i]))
            #       old N[i]s need to be add to list+params
            XlXr[i] = fx[i] * lamda * ( tam1 + tam2 - kW + D * math.cos(omega))
        else:
            XlXr[i] = 0

        try:
            tam3 = (c * beta[i] * math.sin(alpha[i]) + u * beta[i] * math.tan(phi)) / FS
            tam4 = math.cos(alpha[i]) + math.sin(alpha[i]) + math.tan(phi) / FSss
            N[i] = (W[i] + XlXr[i] - tam3) / tam4
        except:
            N[i] = 0
        
        
    return N
     
#--------------------------------------------------------------------------                    

def Cal_FS(Surface, isF, N, c, u, beta, alpha, phi, x, aa, A, R, W):
    tam1 = 0
    tam2 = 0
    tam3 = 0
    for i in range(Surface.shape[0] - 1):
        if isF:
            tam1 = tam1 + c * beta[i] * math.cos(alpha[i])
            tam2 = tam2 + (N[i] - u * beta[i] * math.tan(phi) * math.cos(alpha[i]))
            tam3 = tam3 + N[i] * math.sin(alpha[i]) - A
        else:
            tam1 = tam1 + tam1 + c * beta[i] * R
            tam2 = tam2 + (N[i] - u * beta[i]) * R * math.tan(phi)
            tam3 = tam3 + W[i] * x[i] - A * aa[i]
    print("tam1 = {} ; tam2 = {} ; tam3 = {}".format(tam1, tam2, tam3))
    try:
        FS = (tam1 + tam2) / tam3
        return FS
    except:
        return 0 

#------------------------------------------------------------------------------

def Calculating_FoS(Surface, lamda, fx, center, R, soil, cal, slide, Tolerance):

    # '////////////////////////////////////////////////////
    # '///////////////////////////////////////////////////
    # ' MAIN CALCULATING BELOW
    # ' Building an plot with every surface in vertical of lamda changing series and A list of value F_f,F_m
    # 'After got that, we take intersection point of two lines F_f(lamda) and F_m(lamda).
    # 'At that point, we got F_m = F_f = F.O.S of the current surface with specified center point
    # '   With every lamda values, We use RAPID SOLVER - Multi Iterations to have Convergence Status
    # '   The result of RAPID SOLVER is list of F_m corresponding with lamda and also with F_f
    # '


    # 'MAIN 4 STAGES of FOS CALCULATING
    # 'Main Concept of RAPID SOLVER - Newton-Raphson Algorithm
    # '-----------------------------------------------------------------
    # '   STAGE 1: 1 Iteration : Set X = E = 0 at every slice ; set lamda = 0 ;
    # '   STAGE 2: 4-6 Iterations :    X = E = 0 as XlXr = False in Calculating_N ; to convergence F_m and F_f
    # '   STAGE 3: Main RAPID SOLVER - According to Newton-Raphson Algorithm
    # '       set initial lamda = 2/3 slope between crest and toe
    # '       Selected Tolerance
    # '       Iterations until we got Abs(F_f - F_m) < Tolerance
    # '   STAGE 4: We got a list of F_f(lamda) and F_m(lamda)
    # '       Aferthat, we plot F_f(lamda) and F_m(lamda) ; the intersection is value of F.O.S that we finding
    
    #       DECLARATION AT THE BEGINNING POINT
    #       lamda : np array
    N = np.array([0] * (Surface.shape[0] - 1))
    FS = 1
    FS_f = np.array([float(0)] * lamda.shape[0])
    FS_m = np.array([float(0)] * lamda.shape[0])

    for i in range(lamda.shape[0]):
        #   STAGE  1 : No X and E ; lamda = 0
        FS_f_crr = Cal_FS(Surface, True, N, soil.c, cal.u, slide.beta, slide.alpha, soil.phi, slide.x, slide.a, cal.A, R, slide.W)
        FS_m_crr = Cal_FS(Surface, False,N, soil.c, cal.u, slide.beta, slide.alpha, soil.phi, slide.x, slide.a, cal.A, R, slide.W)
        for j in range(fx.shape[0]):
            #   STAGE 2
            for k in range(6):
                N = Cal_N(Surface, False, FS, lamda[i], fx, soil.c, cal.u, slide.alpha, slide.beta, soil.phi, slide.W, cal.omega, cal.kW)
                FS_f_crr = Cal_FS(Surface, True, N, soil.c, cal.u, slide.beta, slide.alpha, soil.phi, slide.x, slide.a, cal.A, R, slide.W)

            for k in range(6):
                N = Cal_N(Surface, False, FS, lamda[i], fx[j], soil.c, cal.u, slide.alpha, slide.beta, soil.phi, slide.W, cal.omega, cal.kW)
                FS_m_crr = Cal_FS(Surface, False, N, soil.c, cal.u, slide.beta, slide.alpha, soil.phi, slide.x, slide.a, cal.A, R, slide.W)

            #   STAGE 3
            count = 0
            while (abs(FS_f_crr - FS_m_crr) > Tolerance) or (count > 2000):
                count += 1
                N = Cal_N(Surface, True, FS_f_crr, lamda[i], fx[j], soil.c, cal.u, slide.alpha, slide.beta, soil.phi, slide.W, cal.omega, cal.kW)
                FS_f_crr = Cal_FS(Surface, True, N, soil.c, cal.u, slide.beta, slide.alpha, soil.phi, slide.x, slide.a, cal.A, R, slide.W)
                N = Cal_N(Surface, True, FS_m_crr, lamda[i], fx[j], soil.c, cal.u, slide.alpha, slide.beta, soil.phi, slide.W, cal.omega, cal.kW)
                FS_m_crr = Cal_FS(Surface, False, N, soil.c, cal.u, slide.beta, slide.alpha, soil.phi, slide.x, slide.a, cal.A, R, slide.W)
        FS_f[i] = FS_f_crr
        FS_m[i] = FS_m_crr

    #   STAGE 4:
    temp = 0
    for i in range(lamda.shape[0] - 1):
        if (FS_f[i] - FS_f[i + 1]) * (FS_m[i] - FS_m[i + 1]) < 0:
            temp = i
            break
    FS = line_intersection(((FS_f[temp], lamda[i]), (FS_f[temp + 1], lamda[temp + 1])), ((FS_m[temp], lamda[temp]), (FS_m[temp + 1], lamda[temp + 1])))
    
    return FS




#===========================================================================
def Cal_BW(Surface, center, R, dx):
    isLeft, first = index_first(Surface)
    if isLeft:
        xi = center[0] - math.sqrt(R**2 - center[1]**2)
    else:
        xi = center[0] + math.sqrt(R**2 - center[1]**2)
    return abs(xi - first) * dx
 

