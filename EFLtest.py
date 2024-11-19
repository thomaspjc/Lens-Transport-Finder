#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 09:04:34 2024

Testing the EFL of multiple lens systems

@author: thomas
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi as pi
import EFLTools as Tools
from tqdm import tqdm
from ThorLabsLenses import WavelengthAdapter
from FusedSilica import FusedSilica, CalciumFluoride

def FullAnalysis(lenses, d1):
    #--- Determine the true focal lengths of the lenses ---
    trueFocals = LensAdapter(lenses)
    f1, f2 = trueFocals
    
    #--- Extract the specifics of the setup ---
    return EFLTest(f1, f2, d1)
    
    

def EFLTest(f1, f2, d1):
    #Intialising the q parameter
    w0 = 4e-3
    wavelength = 253e-9
    q0 = 1j* pi/wavelength * w0**2
    
    #Fidning the matrices to propagate
    F = propagator1(f1, f2, d1)
    
    #Using the q parameter to determine the distance from f2 to the focal point
    qF2 = (F[0][0]*q0 + F[0][1])/(F[1][0]*q0 + F[1][1])
    d2 = -qF2.real
    
    M = propagator2(f1, f2, d1, d2)
    # --- Finding the Gaussian Width of the beam at the focal point ---
    qf = (M[0][0]*q0 + M[0][1])/(M[1][0]*q0 + M[1][1])
    finalRadius = np.abs((-pi/wavelength * np.imag(1/qf))**(-1/2))
    
    
    # --- Finding the charactersitics of the setup ---
    
    # --- Preparing the Matrix Representations --- 
    A = 1-d1/f1
    B = d1
    C = - (f1+f2-d1)/(f1*f2)
    D = 1-d1/f2
    
    # --- Determine the EFL for a 2 Lens system ---
    EFL = 1 / (1/f1 + 1/f2 - d1/(f1 * f2))

    # --- Determine the Principal Plane distance from the last lens ---
    H2 = (1-A)/C
    
    return d2, finalRadius, EFL, H2

def propagator1 (f1, f2, d1):
    F1 = np.array([[1, 0], [-1/f1, 1]])
    F2 = np.array([[1,0], [-1/f2, 1]])
    D1 = np.array([[1,d1], [0,1]])

    #Determining the full focal matrix & the full system matrix
    F = F2 @ D1 @ F1
    return F
def propagator2 (f1, f2, d1, d2):
    F1 = np.array([[1, 0], [-1/f1, 1]])
    F2 = np.array([[1,0], [-1/f2, 1]])
    D1 = np.array([[1,d1], [0,1]])
    D2 = np.array([[1,d2], [0,1]])
    
    #Determining the full focal matrix & the full system matrix
    F = F2 @ D1 @ F1
    M = D2 @ F
    
    return M

def LensAdapter(lenses):
    """
    

    Parameters
    ----------
    lenses : list
        [Radius of Curvature, Material], [Radius of Curvature, Material]

    Returns
    -------
    effLens1 : TYPE
        DESCRIPTION.
    effLens2 : TYPE
        DESCRIPTION.

    """ 
    trueFocal = []
    for lens in lenses:
        #finding the index of refraction of the material
        n = index(lens[1])
        if lens[1] == "UVFS":
            focal = lensMaker(n, lens[0])
        elif lens[1] == 'CaF2':
            focal = lensMaker(n, lens[0])
        else:
            print('The material is not available')
            return

        trueFocal += [focal]
        
    return trueFocal

def lensMaker(n, R1, R2 = np.inf):
    return 1 / ((n - 1) * (1 / R1 - 1 / R2))

def index(material, wavelength = 253):
    #the wavelength is in nm within the wavelength addapter code
    if material == "UVFS":
        n = FusedSilica(wavelength)
    if material == 'CaF2':
        n = CalciumFluoride(wavelength)
    return n 


if __name__ == "__main__":
    """
    Example Lenses
    -100mm PLCC UVFS ThorLabs: Radius of Curvature -46.0 mm
    175mm PLCX UVFS ThorLabs: Radius of Curvature 80.5 mm
    """
    lenses = [[-46e-3, "UVFS"], [80.5e-3, "UVFS"]]
    separation = 92.33e-3 #distance between the lenses
    d2, finalRadius, EFL, H2 = FullAnalysis(lenses, separation)
    
    print("Distance from the second lens to the cathode:", d2,
          "\nDistance from the second lens to the wheel:", d2-1.200,
          "\nFinal Radius from Gaussian Analysis", finalRadius,
          "\nEFL of the 2 lens system", EFL,
          "\nH2 Plane of the system", H2)







