#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 17:26:05 2024

Transport Lens Finder using Geometric Optics 

@author: thomas
"""
import sympy as sp
import matplotlib.pyplot as plt
import numpy as np

# --- Globals --- 
Lenses = [0.1, 0.150, 0.2]
maxHeight = 0.0508 # 2inches
dConstraints = [0, 10, 0, 30, 0, 1.4] #aMin, aMax, bMin, bMax, cMin, cMax
d0 = 0 #Setting the initial distance from the first lens 
v0 = np.array([4e-3, 4e-7]) 
inputBeam = [[1, d0], [0, 1]]


def propagatorMatrix(z):
    """
    Generates the array to propagate the beam in free space

    Parameters
    ----------
    z : float
        distance to be propagated in meters

    Returns
    -------
    array
        propagation array

    """
    return np.array([[1, z], [0, 1]])


def lensMatrix(f):
    """
    Generates the array to change the angle of the beam at a lens

    Parameters
    ----------
    f : float
        focal length of the lens

    Returns
    -------
    array
        lens array

    """
    return np.array([[1, 0], [-1/f, 1]])


def TransportFinder(lenses, magnification, distantConstraints = dConstraints, maxHeight = 0.0508, inputBeam = inputBeam, inputVector = v0):
    """
    Determines possible setups for a 3 lens system given a specified magnification
    Possible solutions will be within the distant constraints given
    The beam will be less than the specified max height at the lenses

    Parameters
    ----------
    lenses : array
        contains the 3 focal lengths in order of application f1 -> f2 -> f3
    magnification : float
        magnification of the optical system
    distantConstraints : list, optional
        distant constraints as specified in the gloabal, The default is dConstraints.
    maxHeight : float, optional
        The maximum width of the beam at each lens, The default is 0.0508.
    inputBeam : np.array 2x2, optional
        The initial propagation of the beam before the first lens, The default is inputBeam.
    inputVector : 2x1 array, optional
       The input vector of the beam before inputBeam. The default is v0.

    Returns
    -------
    parametric : set of parametric equations 
    
    solutions : solutions to the parametric equations
    
    [A, B, C] : each A[i], B[i], C[i] is a possible setup of lenses 
            Could take the np.transpose to simplify reading

    """
    # --- Extracting Variables --- 
    f1, f2, f3 = lenses
    
    # --- Define symbols for distances --- 
    a, b, c = sp.symbols('a b c')
    
    # --- Creating the matrices for propagation --- 
    Fs = []
    for f in [f1, f2, f3]:
        Fs.append(lensMatrix(f))
    F1, F2, F3 = Fs

    Ds = []
    for z in [a, b, c]:
        Ds.append(propagatorMatrix(z))
    A, B, C = Ds #A, B, C = [[1, a], [0, 1]], [[1, b], [0, 1]], [[1, c], [0, 1]]
    
    # --- Using Geometric Optics to determine the final transfer matrix ---
    transferMatrices = [C, F3, B, F2, A, F1, inputBeam]
    
    transferFinal = np.eye(2)
    for matrix in transferMatrices:
        transferFinal = np.dot(transferFinal, matrix)
    
    # --- Building the equations that we need to solve --- 
    equations = [sp.Eq(transferFinal[0][0], magnification),
                 sp.Eq(transferFinal[0][1], 0),
                 sp.Eq(transferFinal[1][1], 1/magnification)]

    
    solutions = sp.solve(equations, (a, b, c))
    
    # --- Creating parametric function of the solutions ---
    aExpr, bExpr, cExpr = solutions[0] #Extracting the symbolic functions
    #Parametrizing the functions
    aFunc = sp.lambdify(b, aExpr, "numpy")
    bFunc = sp.lambdify(b, bExpr, "numpy")
    cFunc = sp.lambdify(b, cExpr, "numpy")
    parametric = [aFunc, bFunc, cFunc]
    
    # --- Selecting only values that are within the constraints ---
    try : 
        A, B, C = rangeFilter(parametric, distantConstraints, maxHeight, transferMatrices, v0, points = 1e4)
        plt.plot(B, C, '--')
        plt.plot(B, A, '--')
        plt.show()
    except TypeError:
        print("There are no solutions with the given contraints and lenses")
        pass
        return False
    
    return parametric, solutions, [A, B, C]

def rangeFilter(parametric, distantConstraints, maxHeight, transferMatrices, v0, points = 1e4):
    """
    Function that filters the solutions using the distant constraints and the maxHeight parameter

    Parameters
    ----------
    parametric : TYPE
        DESCRIPTION.
    distantConstraints : TYPE
        DESCRIPTION.
    maxHeight : TYPE
        DESCRIPTION.
    transferMatrices : TYPE
        DESCRIPTION.
    v0 : TYPE
        DESCRIPTION.
    points : TYPE, optional
        DESCRIPTION. The default is 1e4.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # --- Extracting Parameters --- 
    aFunc, _ , cFunc = parametric
    aMin, aMax, bMin, bMax, cMin, cMax = distantConstraints
    bSpace = np.arange(bMin, bMax+1, (bMax+1 - bMin)/1e4)
    _, F3, _, F2, _, F1, inputBeam = [np.array(matrix) for matrix in transferMatrices]
    v0 = np.array(v0)
    
    # --- Filtering based on distance constraints ---
    #Defining values for the a and c distances from the b space
    A = aFunc(bSpace)
    B = bSpace
    C = cFunc(bSpace)
    
    # Create a boolean mask based on the constraints
    mask = ((A >= aMin) & (A <= aMax) &
            (C >= cMin) & (C <= cMax))
    #Filtering the values using the boolean mask 
    A = A[mask]
    B = B[mask]
    C = C[mask]
    
    # --- Filtering based on the beam waist at the lenses in each setup --- 
    #initialising the new mask for efficiency
    heightMask = np.ones_like(A, dtype = bool)
    
    
    for i, (a,b,c) in enumerate(zip(A, B, C)):
        #Creating the transport matrices for each setup
        D1 = propagatorMatrix(a)
        D2 = propagatorMatrix(b)
        D3 = propagatorMatrix(c)
        
        vectorA = F2 @ D1 @ F1 @ inputBeam @ v0
        vectorB = F3 @ D2 @ vectorA
        vectorC = D3 @ vectorB #Could remove this since we have a defined M
    
        heightA, heightB, heightC = vectorA[0], vectorB[0], vectorC[0]
        
        if (abs(heightA) > maxHeight or abs(heightB) > maxHeight or abs(heightC > maxHeight)):
            heightMask[i] = False
        
    #Applying the beam waist filter to the possible setups
    A = A[heightMask]
    B = B[heightMask]
    C = C[heightMask]
    
    if len(A) == 0:
        return False
    
    return A, B, C


def waistPlotting(solution, lenses, inputBeam = inputBeam, inputVector = v0):
    """
    

    Parameters
    ----------
    solution : TYPE
        DESCRIPTION.
    lenses : TYPE
        DESCRIPTION.
    inputBeam : TYPE, optional
        DESCRIPTION. The default is inputBeam.
    inputVector : TYPE, optional
        DESCRIPTION. The default is v0.

    Returns
    -------
    None.

    """
    # --- Extracting Parameters ---
    A, B, C = solution
    d0 = inputBeam[0][1]
    f1, f2, f3 = lenses
    # --- Setting up the plots ---
    #Propagation space between each lens
    z0 = np.linspace(0, d0, int(1e1))
    z1 = np.linspace(0, A, int(1e1))
    z2 = np.linspace(0, B, int(1e1))
    z3 = np.linspace(0, C, int(1e1))
    
    # --- Preparing the matrices --- 
    Fs = [lensMatrix(f) for f in [f1, f2, f3]]
    F1, F2, F3 = Fs
    
    M0 = [propagatorMatrix(z) for z in z0]
    M1 = [propagatorMatrix(z) for z in z1]
    M2 = [propagatorMatrix(z) for z in z2]
    M3 = [propagatorMatrix(z) for z in z3]

    #Assume that each seperation has the same number of points
    vI = [np.dot(F1 @ M0[i], v0) for i in range(len(M0))]
    va = [np.dot(F2 @ M1[i], vI[-1]) for i in range(len(M1))]
    vb = [np.dot(F3 @ M2[i], va[-1]) for i in range(len(M2))]
    vc = [M3[i] @ vb[-1] for i in range(len(M3))]
    
    # --- Extract the beam height from each vector --- 
    h0 = [vI[i][0] for i in range(len(vI))]
    ha = [va[i][0] for i in range(len(va))]
    hb = [vb[i][0] for i in range(len(vb))]
    hc = [vc[i][0] for i in range(len(vc))]
    
    
    # --- Plotting ---
    plt.plot(z0, h0, 'b-')
    plt.plot(z1+d0, ha, 'b-')
    plt.plot(z2+d0+A, hb, 'b-')
    plt.plot(z3+d0+A+B, hc, 'b-')
    plt.plot(z0, (-1 * np.array(h0)), 'b-')
    plt.plot(z1+d0, (-1 * np.array(ha)), 'b-')    
    plt.plot(z2+d0+A, (-1 * np.array(hb)), 'b-')
    plt.plot(z3+d0+A+B, (-1 * np.array(hc)), 'b-')
    
    # Function to draw a lens
    def drawLens(position, length, label = False):
        # Draw the vertical line
        if label:
            plt.plot([position, position], [-length, length], 'k--', linewidth=2, label = 'Lens')
        else:
            plt.plot([position, position], [-length, length], 'k--', linewidth=2)
        
    drawLens(d0, maxHeight, label = True)
    drawLens(A+d0, maxHeight)
    drawLens(B+A+d0, maxHeight)

    plt.legend()
    plt.title("Beam Waist through the Optical Setup")
    plt.xlabel("Distance (m)")
    plt.show()
    return
    
    
"""Ms = np.arange(0.5, 1.5, 0.5)
for i in Ms:
    solutions = TransportFinder(Lenses, i)
    if solutions != False:
        print(i)
"""
    
solutions = TransportFinder(Lenses, 1)

selected = np.transpose(solutions[-1])[150]
waistPlotting(selected, Lenses)




































