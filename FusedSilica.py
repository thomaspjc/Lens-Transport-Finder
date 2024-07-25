#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 16:49:55 2024

@author: thomas
"""

import numpy as np
import matplotlib.pyplot as plt

def calculate_refractive_index(wavelength_nm):
    # Convert wavelength from nm to micrometers
    wavelength_um = wavelength_nm / 1000.0
    lambda2 = wavelength_um**2
    term1 = 0.6961663 * lambda2 / (lambda2 - 0.0684043**2)
    term2 = 0.4079426 * lambda2 / (lambda2 - 0.1162414**2)
    term3 = 0.8974794 * lambda2 / (lambda2 - 9.896161**2)
    n2_minus_1 = term1 + term2 + term3
    n = np.sqrt(n2_minus_1 + 1)
    return n

def plot_refractive_index():
    # Define wavelength ranges
    wavelengths_full = np.linspace(200, 7000, 1000)
    wavelengths_zoomed = np.linspace(200, 800, 1000)
    
    # Calculate refractive indices
    refractive_indices_full = calculate_refractive_index(wavelengths_full)
    refractive_indices_zoomed = calculate_refractive_index(wavelengths_zoomed)
    
    # Calculate the index of refraction for 253 nm
    wavelength_253 = 253
    n_253 = calculate_refractive_index(wavelength_253)
    wavelength_ThorLabs = 588
    nThorLabs = calculate_refractive_index(wavelength_ThorLabs)
    
    # Create the main plot
    fig, ax_main = plt.subplots(figsize=(12, 8))
    
    # Plot the zoomed range (200 - 800 nm) in the main plot
    ax_main.plot(wavelengths_zoomed, refractive_indices_zoomed, label='200 - 800 nm', color='red')
    
    ax_main.text(255, 1.456, f'n at 253nm = {n_253:.6f}\nn at 588nm = {nThorLabs:.6f}', fontsize=12, ha='center', va='bottom',
                 bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    
    # Labels and title for the main plot
    ax_main.set_xlabel('Wavelength (nm)')
    ax_main.set_ylabel('Refractive Index (n)')
    ax_main.set_title('Refractive Index of Fused Silica in Terms of Wavelength (nm)')
    ax_main.legend(loc = 'lower left')
    ax_main.grid(True)
    
    # Create the inset plot for the full range (200 - 6000 nm)
    ax_inset = fig.add_axes([0.59, 0.545, 0.3, 0.3])  # [left, bottom, width, height]
    
    # Plot the full range in the inset
    ax_inset.plot(wavelengths_full, refractive_indices_full, label='200 - 6000 nm', color='blue')
    
    # Labels and title for the inset plot
    ax_inset.set_xlabel('Wavelength (nm)')
    ax_inset.set_ylabel('Refractive Index (n)')
    ax_inset.set_title('Full Range')
    ax_inset.grid(True)
    
    # Show the plot
    plt.show()

# Example usage
#plot_refractive_index()

