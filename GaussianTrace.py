#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 09:18:06 2024

2D Beam Propagation Plotting

@author: thomas
"""

import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi as pi
from tqdm import tqdm
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
from cmocean import cm
import seaborn as sns
from matplotlib.colors import Normalize


# --- Globals ---
wavelength = 253 * 1e-9
w0 = 4
f = 1.2
k0 = 2* pi / wavelength
#extent = [-8 * w0, 8 * w0]
z0 = pi/wavelength * w0**2
extent = [-1.27 * 1e-2, 1.27 * 1e-2]

def fluka_color_map(self):
        
        flukacolors = [(1.0, 1.0, 1.0), (0.9, 0.6, 0.9), (1.0, 0.4, 1.0), (0.9, 0.0, 1.0), (0.7, 0.0, 1.0), (0.5, 0.0, 0.8), (0.0, 0.0, 0.8),
                   (0.0, 0.0, 1.0), (0.0, 0.6, 1.0), (0.0, 0.8, 1.0), (0.0, 0.7, 0.5), (0.0, 0.9, 0.2), (0.5, 1.0, 0.0), (0.8, 1.0, 0.0),
                   (1.0, 1.0, 0.0), (1.0, 0.8, 0.0), (1.0, 0.5, 0.0), (1.0, 0.0, 0.0), (0.8, 0.0, 0.0), (0.6, 0.0, 0.0), (0.0, 0.0, 0.0)]
        cmap_name = 'fluka'

        return LinearSegmentedColormap.from_list(cmap_name, flukacolors, N=300)


def Transfer(inputQ, transferMatrix):
    
    # --- Extracting the Parameters ---
    (A,B),(C,D) = transferMatrix[0], transferMatrix[1]
    # --- Applying the transfer Matrix ---
    qPrime = (A * inputQ + B)/(C*inputQ + D)
    
    # --- Extracting new parameters for the beam ---
    w_z = np.abs((-pi/wavelength * np.imag(1/qPrime))**(-1/2))
    R_z = ((1/qPrime).real)**(-1)
    
    return qPrime, w_z, R_z
    
def Lens(f):
    return np.array([[1,0],[-1/f, 1]])

def Propagate(z):
    return np.array([[1, z],[0, 1]])


def BeamPlotting(distances, focals, inputBeam, sampling):
    
    BuildingInput() #This function finds the inputQ for a given width and divergence
    Qs, Ws = waistIdentification()
    intensityMatching()
    return



def waistIdentification(distances, focals, inputQ, sampling):
    
    # --- Extracting the data ---
    A, B, C = (Propagate(d) for d in distances)
    lenses = [Lens(f) for f in focals]
    
    totalTravel = 0
    totalTransfer = np.eye(2)
    Qs = [[inputQ]]
    Ws = [w0]
    for i, d in tqdm(enumerate(distances), desc = 'Iterating Through Optical System'):
        
        Z = np.arange(0, d + sampling, sampling)
        D = [Propagate(z) for z in Z]
       
        #Create an array of transfer matrices to apply to 
        transferMatrix = [d @ lenses[i] @ totalTransfer for d in D]

        #Extracting the parameters for each iteration of the propagated gaussian
        transfered = [Transfer(Qs[0][0], T) for T in transferMatrix]
        qs, ws, _ = np.array(transfered).T
        #Preparing for next iteration
        totalTravel += d
        totalTransfer = D[-1] @ lenses[i] @ totalTransfer 
        
        Qs += [qs]
        Ws += [ws.real]
    return Qs[1:], Ws[1:]


def intensityMatching(Qs, inputBeam, k0 = k0, ySampling = 1e4):
    #Make the intensity pattern for a 1D Gaussian using the Q parameter for every sampling step
    qArray = np.concatenate(Qs[:])

    inputWidth, _ = inputBeam
    ySpace = np.linspace(-8*inputWidth, 8*inputWidth, int(ySampling))
    
    qq, yy = np.meshgrid(qArray, ySpace, indexing= 'xy')
    gaussianCut = 1/np.sqrt(qq) * np.exp(-1j * k0 * yy**2 / (2 * qq))
    
    return gaussianCut

def BuildingInput(inputBeam, wavelength = 253e-9):
    wZ, alpha = inputBeam
    z = wZ / np.tan(2*alpha)
    w0 = np.sqrt(wZ**2 + np.sqrt(pi**2 * wZ**4 - 4 * z**2 * wavelength**2) / pi) / np.sqrt(2)
    w0 = np.sqrt(wZ**2 - np.sqrt(pi**2 * wZ**4 - 4 * z**2 * wavelength**2) / pi) / np.sqrt(2)
    
    z0 = pi * w0**2/wavelength
    qZ = z + 1j * z0
    
    return qZ





# --- Initialization ---
v0 = np.array([4, 4e-4])
distances = [114.7, 10.7, 1200]
focals = np.array([-91.52910582065768, -18210.995922289923, 181.4767152124487])
inputBeam = [4, 4e-4]
#inputQ = -200 + 1j * pi / wavelength * w0**2
inputQ = BuildingInput(inputBeam)
sampling = 0.1
ySampling = 1e4



# --- Waist Identification ---
Qs, Ws = waistIdentification(distances, focals, inputQ, sampling)   
print(list(Qs[0][0:5]))
# --- Concatenate and Plot Beam Widths ---
ws = np.concatenate(Ws)

# --- Intensity Matching and Plot ---
norm = Normalize(vmin=0, vmax=25 / ySampling)
gaussians = intensityMatching(Qs[:], inputBeam, ySampling = ySampling)
data = np.abs(gaussians)**2
totalIntensity = np.sum(data, axis=0)
cmaps = [cm.curl, cm.thermal, cm.dense, 'mako', cm.solar, cm.haline, 'GnBu_r']
plt.imshow(data / totalIntensity, aspect='auto', cmap=cmaps[0],
           extent=[0, np.sum(distances), -8 * inputBeam[0], 8 * inputBeam[0]], norm = norm)
plt.colorbar()

# --- Plot the Widths on Top of the Intensity Plot ---

z_positions = np.linspace(0, np.sum(distances), len(ws))
colors = ['white', 'red', 'white']
start_idx = 0

for i, color in enumerate(colors):
    end_idx = start_idx + len(Ws[i])
    plt.plot(z_positions[start_idx:end_idx], Ws[i], color=color, linestyle='--',
             linewidth=1.3, label=r'Width ($1/e^2$)' if i == 0 else "")
    plt.plot(z_positions[start_idx:end_idx], -Ws[i], color=color, linestyle='--', linewidth=1.3)
    start_idx = end_idx
    
# --- Adding a text box ---
textstr = f'Final Beam Height = {Ws[-1][-1]:.6g}mm'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(400, -28, textstr, fontsize=10,
        verticalalignment='top', bbox=props)

# --- Add Labels, Title, and Legend ---
plt.xlabel('Propagation distance (mm)')
plt.ylabel('Transverse position (mm)')
plt.title('Beam Propagation through 3 Lens System')
plt.legend()
plt.show()
    
