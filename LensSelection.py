#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 10:24:01 2024

Lens Selection

@author: thomas
"""
import numpy as np
from tqdm import tqdm

# --- Globals ---
fConcave = [-25, -50, -75, -100] #Fused Silica 1 inch in mm #uncoated UV
fConvex1Inch = [35, 40, 50, 75, 100, 125, 150, 200, 250, 300, 500, 750, 1000]# Fused Silica 1 inch in mm uncoated UV
fConvex2Inch = [60,75,100,150,200,250,300,500,1000]# Fused Silica 2 inch in mm uncoated UV


def select2Lens(F1, F2, mConstraint, dConstraint):
    # Create all possible pairs
    f1, f2 = np.meshgrid(fConcave, fConvex1Inch, indexing='ij')
    f1 = f1.ravel()
    f2 = f2.ravel()
    
    # Calculate conditions
    sum_condition = f1 + f2 < 100
    magnification = -f2 / f1
    magnification_condition = (magnification > 2) & (magnification < 2.6)
    
    # Combine conditions
    valid_combinations = sum_condition & magnification_condition
    #print(valid_combinations)
    #print(len(f1), len(f2))
    # Extract valid pairs
    valid_f1 = f1[valid_combinations]
    valid_f2 = f2[valid_combinations]
    valid_pairs = list(zip(valid_f1, valid_f2))
    
    # Print the valid pairs
    print("Valid pairs (f1, f2):")
    for pair in valid_pairs:
        print(pair)
    return valid_pairs
    

    
    
def Effective3Lens(F1, F3, d1, d2, inputBeam, fEff = 1200, thetaPlate  = 0.057):
    """
    Finds the arangement of lenses and their positions such that a phase plate 
    is the second of a 3 lens set up and we find their effective focal length 

    Parameters
    ----------
    F1 : TYPE
        DESCRIPTION.
    F3 : TYPE
        DESCRIPTION.
    d1 : TYPE
        DESCRIPTION.
    d2 : TYPE
        DESCRIPTION.
    inputBeam : TYPE
        DESCRIPTION.
    fEff : TYPE, optional
        DESCRIPTION. The default is 1200.
    thetaPlate : TYPE, optional
        DESCRIPTION. The default is 0.057.

    Returns
    -------
    valid_pairs : TYPE
        DESCRIPTION.

    """
    inch1 = 25.4 
    inch2 = 50.8 
    # --- Extracting variables --- 
    x, alpha = inputBeam
    
    
    #Creating a 4D Meshgrid to find all possible combinations of lens arrangements
    f1, f3, d1, d2 = np.meshgrid(F1, F3, d1, d2, indexing ='xy')

    # --- Determine the focal length of the phase plate ---
    #Determine the Diameter of the beam at the phase plate
    f2Diameter = (x - d1 * x / f1 + d1 * alpha) *2
    #Determine the absolute value of the focal length 
    f2abs = f2Diameter / (2 * np.tan(thetaPlate/2 * np.pi/180))
    f2 = - f2abs

    
    # --- Determine the Diameter at the final lens --- 
    f3Radius = ((f1 * f2 * x + d1 * (d2 - f2) * (x - f1 * alpha) -
                  d2 * (f1 * x + f2 * x - f1 * f2 * alpha)) / (f1 * f2))
    

    
    # --- Creating the mask to select functional arrangements --- 
    mask =(
        (fEff == np.round(((f1 * f2 * f3) / (f1 * f2 - d2 * (f1 + f2) +
                                   d1 * (d2 - f3) + f1 * f3 + f2 * f3)), -1))
        & 
        (f2Diameter < 0.95 * inch1) & (f2Diameter > 0.8 * inch1)
        &
        (f3Radius < inch2/2)
        &
        (d1+d2 < 500)
        )
    d3 = (np.round(((f1 * f2 * f3) / (f1 * f2 - d2 * (f1 + f2) +
                               d1 * (d2 - f3) + f1 * f3 + f2 * f3)), 0))
    #print(d3)
    # --- Selecting items which correspond to valid solution pairs --- 
    #Updating the progress bar

    f1Valid = f1[mask]
    f2Valid = f2[mask]
    f3Valid = f3[mask]
    d1Valid = d1[mask]
    d2Valid = d2[mask]
    d3Valid = d3[mask]
    #Optionals which are recovered through ray tracing
    f2Diameter = f2Diameter[mask]
    f3Radius = f3Radius[mask]
    
    
    # --- Pairing up the pair wise valid solutions --- 
    valid_pairs = list(zip(f1Valid, f2Valid, f3Valid, d1Valid, d2Valid, d3Valid))
    return valid_pairs

def Finder3Lens(F1, F3, d1, d2, inputBeam, d3 = 1200, thetaPlate  = 0.057,
                tolerance = 1e-1):
    """
    Finds the arangement of lenses and their positions such that a phase plate 
    is the second of a 3 lens set up and we find their effective focal length 

    Parameters
    ----------
    F1 : TYPE
        DESCRIPTION.
    F3 : TYPE
        DESCRIPTION.
    d1 : TYPE
        DESCRIPTION.
    d2 : TYPE
        DESCRIPTION.
    inputBeam : TYPE
        DESCRIPTION.
    fEff : TYPE, optional
        DESCRIPTION. The default is 1200.
    thetaPlate : TYPE, optional
        DESCRIPTION. The default is 0.057.

    Returns
    -------
    valid_pairs : TYPE
        DESCRIPTION.

    """
    inch1 = 25.4 
    inch2 = 50.8 
    # --- Extracting variables --- 
    x, alpha = inputBeam
    
    
    #Creating a 4D Meshgrid to find all possible combinations of lens arrangements
    f1, f3, d1, d2 = np.meshgrid(F1, F3, d1, d2, indexing ='xy')

    # --- Determine the focal length of the phase plate ---
    #Determine the Diameter of the beam at the phase plate
    f2Diameter = (x - d1 * x / f1 + d1 * alpha) *2
    #Determine the absolute value of the focal length 
    f2abs = f2Diameter / (2 * np.tan(thetaPlate/2 * np.pi/180))
    f2 = -f2abs

    
    # --- Determine the Diameter at the final lens --- 
    f3Radius = ((f1 * f2 * x + d1 * (d2 - f2) * (x - f1 * alpha) -
                  d2 * (f1 * x + f2 * x - f1 * f2 * alpha)) / (f1 * f2))
    
    finalRadius = ((1 / (f1 * f2 * f3)) * (-d3 * f1 * f2 * x - d3 * f1 * f3 *
                    x - d3 * f2 * f3 * x + f1 * f2 * f3 * x + d3 * f1 * f2 *
                    f3 * alpha + d1 * (-f2 * f3 + d2 * (-d3 + f3) + d3 * 
                    (f2 + f3)) * (x - f1 * alpha) + d2 * (d3 - f3) * 
                    (f2 * x + f1 * (x - f2 * alpha)))
                   )

    # --- Creating the mask to select functional arrangements --- 
    mask =(
        (tolerance > np.abs(finalRadius))
        & 
        (f2Diameter < 0.82 * inch1) & (f2Diameter > 0.71 * inch1)
        &
        (f3Radius < 0.82 * inch1/2)
        &
        (d1+d2 < 200)
        )
    d3 = (np.round(((f1 * f2 * f3) / (f1 * f2 - d2 * (f1 + f2) +
                               d1 * (d2 - f3) + f1 * f3 + f2 * f3)), 0))
    #print(d3)
    # --- Selecting items which correspond to valid solution pairs --- 
    #Updating the progress bar

    f1Valid = f1[mask]
    f2Valid = f2[mask]
    f3Valid = f3[mask]
    d1Valid = d1[mask]
    d2Valid = d2[mask]
    d3Valid = d3[mask]
    #Optionals which are recovered through ray tracing
    f2Diameter = f2Diameter[mask]
    f3Radius = f3Radius[mask]
    finalValidRadius = finalRadius[mask]
    
    # --- Pairing up the pair wise valid solutions --- 
    valid_pairs = list(zip(f1Valid, f2Valid, f3Valid, d1Valid, d2Valid, finalValidRadius, f2Diameter/2, f3Radius))
    return valid_pairs



























