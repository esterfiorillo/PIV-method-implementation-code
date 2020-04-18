# -*- coding: utf-8 -*-
"""

 (               )     )   |
 )\ )  *   )  ( /(  ( /(   |
(()/(` )  /(  )\()) )\())  |
 /(_))( )(_))((_)\ ((_)\   | LTHN PIV: In an opensource PIV toolbox
(_)) (_(_())  _((_) _((_)  |
| |  |_   _| | || || \| |  |
| |__  | |   | __ || .` |  | Website: https://github.com/esterfiorillo/PIV-method-implementation-code
|(___|(|_|   |_||_||_|\_|  |
                           | 
 )\ ) )\ )                 | CDTN - Centro de Desenvolvimento da Tecnologia Nuclear
(()/((()/( (   (           | LTHN - Laboratório de Termo-Hidráulica e Neutrônica
 /(_))/(_)))\  )\          | Belo Horizonte, MG, Brasil
(_)) (_)) ((_)((_)         |
| _ \|_ _|\ \ / /          | @authors: esterfiorillo, acampagnole 
|  _/ | |  \ V /           |
|_|  |___|  \_/            |


"""

# License
#     This file is part of LTHN PIV.
#
#     LTHN PIV is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     LTHN PIV is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with LTHN PIV.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     ./processing/interrogation_window.py
#
# Description
#     Class that represents the interrogation windows and the functions applyed to then.
#     Atributes:
#         window1: 2d np_array
#         window2: 2d np_array
#         cross_corr: 2d np_array
#     Functions:
#         normxcorr2()
#         find_peak()
#         gauss_subpixel_peak_position
#
#------------------------------------------------------------------------------

"""
   load basic modules
"""
import sys


"""
   try to load fundamental modules
"""
try:
    import numpy as np
    from numpy import log
    from scipy.signal import fftconvolve
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit()    
 


class interrogation_window:

    
    def __init__(self, w1, w2): #Constructor
        self.window1 = w1
        self.window2 = w2

    
    def normxcorr2(self, mode="full"):
    ########################################################################################
    # Author: Ujash Joshi, University of Toronto, 2017                                     #
    # Based on Octave implementation by: Benjamin Eltzner, 2014 <b.eltzner@gmx.de>         #
    # Octave/Matlab normxcorr2 implementation in python 3.5                                #
    # Details:                                                                             #
    # Normalized cross-correlation. Similiar results upto 3 significant digits.            #
    # https://github.com/Sabrewarrior/normxcorr2-python/master/norxcorr2.py                #
    # http://lordsabre.blogspot.ca/2017/09/matlab-normxcorr2-implemented-in-python.html    #
    ########################################################################################  
    #    Input arrays should be floating point numbers.
    #    :param template: N-D array, of template or filter you are using for cross-correlation.
    #    Must be less or equal dimensions to image.
    #    Length of each dimension must be less than length of image.
    #    :param image: N-D array
    #    :param mode: Options, "full", "valid", "same"
    #    full (Default): The output of fftconvolve is the full discrete linear convolution of the inputs.
    #    Output size will be image size + 1/2 template size in each dimension.
    #    valid: The output consists only of those elements that do not rely on the zero-padding.
    #    same: The output is the same size as image, centered with respect to the ‘full’ output.
    #    :return: N-D array of same dimensions as image. Size depends on mode parameter.
   
    # If this happens, it is probably a mistake
    
        if np.ndim(self.window1) > np.ndim(self.window2) or \
                len([i for i in range(np.ndim(self.window1)) if self.window1.shape[i] > self.window2.shape[i]]) > 0:
                print("normxcorr2: self.window1 larger than IMG. Arguments may be swapped.")

        self.window1 = self.window1 - np.mean(self.window1)
        self.window2 = self.window2 - np.mean(self.window2)

        a1 = np.ones(self.window1.shape)
        
        # Faster to flip up down and left right then use fftconvolve instead of scipy's correlate
        
        ar = np.flipud(np.fliplr(self.window1))
        out = fftconvolve(self.window2, ar.conj(), mode=mode)
   
        self.window2 = fftconvolve(np.square(self.window2), a1, mode=mode) - \
                np.square(fftconvolve(self.window2, a1, mode=mode)) / (np.prod(self.window1.shape))

        # Remove small machine precision errors after subtraction
        
        self.window2[np.where(self.window2 < 0)] = 0

        self.window1 = np.sum(np.square(self.window1))
        out = out / np.sqrt(self.window2 * self.window1)

        # Remove any divisions by 0 or very close to 0
        out[np.where(np.logical_not(np.isfinite(out)))] = 0
   
        self.cross_corr = out

        
    def find_peak (self):    
        # Function that returns the coordinates in x and y of the biggest point (peak)
        # in the 2d array resulted from the cross correlation.
        
        ind = self.cross_corr.argmax()
        s = self.cross_corr.shape[1]
        xx_peak = ind // s
        yy_peak = ind % s
        return xx_peak, yy_peak


    def gauss_subpixel_peak_position (self, x_pico, y_pico):
        # Function that returns the gauss interpolation coordinates
        #
        # In:
        #    x_pico: int
        #         x-coordinate of the peak
        #    y_pico: int
        #         y-coordinate of the peak
        # Out:
        #    subp_peak_position[0]: int
        #        x-coordinate in the cross-correlation matrix resulting from gauss interpolation
        #    subp_peak_position[1]: int
        #        y-coordinate in the cross-correlation matrix resulting from gauss interpolation
        
        if x_pico + 1 >= len(self.cross_corr):
            subp_peak_position = x_pico, y_pico
        elif y_pico +1 >= len(self.cross_corr[0]):
            subp_peak_position = x_pico, y_pico
        else:
           
            c = self.cross_corr[x_pico, y_pico]
            cl = self.cross_corr[x_pico - 1, y_pico]
            cr = self.cross_corr[x_pico + 1, y_pico]
            cd = self.cross_corr[x_pico, y_pico - 1]
            cu = self.cross_corr[x_pico, y_pico + 1]
            if np.any(np.array([c, cl, cr, cd, cu]) < 0):
                subp_peak_position = (((x_pico - 1) * cl + x_pico * c + (x_pico + 1) * cr) / (cl + c + cr),
((y_pico - 1) * cd + y_pico * c + (y_pico + 1) * cu) / (cd + c + cu))
            else:    
                subp_peak_position = (x_pico + ((log(cl) - log(cr)) / (2 * log(cl) - 4 * log(c) + 2 * log(cr))),
y_pico + ((log(cd) - log(cu)) / (2 * log(cd) - 4 * log(c) + 2 * log(cu))))
        return subp_peak_position[0], subp_peak_position[1]        
