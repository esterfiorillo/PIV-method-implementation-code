# -*- coding: utf-8 -*-
"""

 (               )     )   |
 )\ )  *   )  ( /(  ( /(   |
(()/(` )  /(  )\()) )\())  |
 /(_))( )(_))((_)\ ((_)\   | LTHN PIV: Is an opensource PIV toolbox
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
#     ./processing/normal_method.py
#
# Description
#     Method 1
#
#------------------------------------------------------------------------------

from pre_processing.image_par import image_par
from processing.interrogation_window import interrogation_window
import sys
from processing.displacement_map import displacement_map

try:
    import numpy as np
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit() 




class normal_method(image_par):

    
    def __init__ (self, im1, im2, w_size, ovl): #Constructor
        image_par.__init__(self, im1, im2)
        self.x_size = w_size
        self.y_size = w_size
        self.overlap = ovl


    def result_dimensions (self):  
    #Function used to find the dimensions of the result of applying cross correlation in the interrogation windows, according to the interrogation window size and overlap
    #In:
    #    x_im: int
    #        x-dimension of the image
    #    y_im: int
    #        y-dimension of the image
    #    window_size: int
    #        size of the interrogation window that will be used to go through the images
    #    overlap: int
    #        
    #Out:
    #    x_dimension: int
    #        x-dimension of matrix resulting from applying cross-correlation between interrogation windows
    #    y_dimension: int
    #       y-dimension of matrix resulting from applying cross-correlation between interrogation windows
    
        x_dimension = (self.tam_x - self.x_size) // (self.x_size - self.overlap) + 1
        y_dimension = (self.tam_y - self.y_size) // (self.y_size - self.overlap) + 1
        return x_dimension, y_dimension


    def first_iteration (self):
    #Function in which the pair of images is traversed with the interrogation windows while applying the cross correlation between them.
    #In this case, the windows will be shifted normally, with no method to make the process more efficient.
    #This step is required to find the first displacement map that will be used in Multigrid.
    #In:
    #    im1: 2d array
    #        first image
    #    im2: 2d array
    #        second image
    #    tam_x: int
    #        x-dimension of the images
    #    tam_y: int
    #        y-dimension of the images
    #    x_size: int
    #        interrogation window size in x axis
    #    y_size: int
    #        interrogation window size in y axis
    #    overlap: int
    #Out:
    #    dpx1: 2d array
    #        displacement in x
    #    dpy1: 2d array
    #        displacement in y    
    
        size_r_x, size_r_y = self.result_dimensions()
        #create matrix for substitute divisions of im1
        a = np.zeros ((self.x_size, self.y_size))
        #create matrix for substitute search window in im2
        b = np.zeros ((self.x_size, self.y_size))

        dpx1 = np.zeros ((size_r_x, size_r_y))
        dpy1 = np.zeros ((size_r_x, size_r_y))
       
        for i in range (0, size_r_x):
            for j in range (0, size_r_y):
                x1_min = i*(self.x_size - self.overlap)
                x1_max = x1_min + self.x_size

                y1_min = j*(self.y_size - self.overlap)
                y1_max = y1_min + self.y_size

                a = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
                b = np.copy (self.im2[x1_min:x1_max, y1_min:y1_max])
               
                janela = interrogation_window (a , b)
                janela.normxcorr2 ('same')
                
                xx, yy = janela.find_peak()
                xx, yy = janela.gauss_subpixel_peak_position (xx, yy)
                xx -= (self.x_size + self.x_size - 1)//2
                yy -= (self.y_size + self.y_size - 1)//2

                dpx1[i][j], dpy1[i][j] = -xx, yy

        res = displacement_map(dpx1, dpy1)
        return res
