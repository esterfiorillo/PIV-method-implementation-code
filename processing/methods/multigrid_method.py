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
#     ./processing/multigrid_method.py
#
# Description
#     Class that represents the displacement map resulting from the methods.
#     It also have the functions applied to the map for the data validation.
#     Atributes:
#         dp_x: 2d np.narray
#         dp_y: 2d np.narray
#     Functions:
#         neighborhood_median(int ii, int jj)
#         calculate_rm (int ii, int jj)
#         calculate_r (int ii, int jj)
#         neighborhood_median(int ii, int jj)
#         replacement3()
#         resize_dp
#         index_guard()
#
#------------------------------------------------------------------------------

from pre_processing.image_par import image_par
from processing.interrogation_window import interrogation_window
import sys

try:
    import numpy as np
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit() 


class multigrid_method (image_par):
# Method2
    def __init__ (self, im1, im2, w_size, ovl, n_passadas): #Constructor
        image_par.__init__(self, im1, im2)
        self.x_size = w_size
        self.y_size = w_size
        self.overlap = ovl
        self.num_passadas = n_passadas

    def result_dimensions (self):
    #Function used to find the dimensions of the result of applying cross correlation in the interrogation windows, according to the interrogation window size and overlap.
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
   
    def multigrid_method1 (self, dp_anterior):
    #Function that implements the Multigrid method, which is an iterative procedure for traversing images with interrogation windows while applying a cross correlation.
    #It consists of moving windows according to the displacement map found in the previous iteration.
    #At the same time, the size of the interrogation windows is reduced between one iteration and another.
    #In addition, a data validation criterion is applied to the displacement maps.
    #In:
    #    num_passadas: int
    #        number of iterations
    #    x_size: int
    #        interrogation window size in x axis
    #    y_size: int
    #        interrogation window size in y axis
    #    overlap: int
    #    tam_x: int
    #        x-dimension of the images
    #    tam_y: int
    #        y-dimension of the images
    #    dpx_anterior: 2d array
    #        previous displacement in x axis found in the function fisrt_iteration
    #    dpy_anterior: 2d array
    #        previous displacement in y axis found in the function fisrt_iteration
    #    im1: 2d array
    #        first image
    #    im2: 2d array
    #        second image
    #Out:
    #    dpx_anterior: 2d array
    #        displacement in x
    #    dpy_anterior: 2d array
    #        displacement in y
    
    #Reference:
    #F. Scarano & M. L. Riethmuller, Iterative multigrid approach in PIV image processing with discrete window offset, Experiments in Fluids 26 (1999)
        aux = 0
        while (aux != self.num_passadas):
            self.x_size = self.x_size//2
            self.y_size = self.y_size//2
            self.overlap = self.overlap//2
            size_r_x, size_r_y = self.result_dimensions()
            ind1_x, ind1_y, ind2_x, ind2_y = dp_anterior.index_guard()
   
            dp_anterior.resize_dp ()

            dpx = np.zeros((size_r_x, size_r_y))
            dpy = np.zeros((size_r_x, size_r_y))
            c = np.zeros((size_r_x, size_r_y))
            d = np.zeros((size_r_x, size_r_y))
            cont1 = 0
            cont2 = 0
            for i in range (0, size_r_x):
                for j in range (0, size_r_y):
                    x1_min = i* (self.x_size - self.overlap)
                    y1_min = j* (self.x_size - self.overlap)
                    if i == ind2_x[cont1]:
                        if j == ind2_y[cont2]:
                            teste1_indice1 = int (ind1_x[cont1])
                            teste1_indice2 = int (ind1_y[cont2])
   
                            if (np.isnan(dp_anterior.dpx[teste1_indice1][teste1_indice2])):
                                auxliar1 = 0
                            else:
                                auxliar1 = dp_anterior.dpx[teste1_indice1][teste1_indice2]
                            if (np.isnan(dp_anterior.dpy[teste1_indice1][teste1_indice2])):
                                auxliar2 = 0
                            else:
                                auxliar2 = dp_anterior.dpy[teste1_indice1][teste1_indice2]
                            x2_min = int (i*(self.x_size - self.overlap) + auxliar1)
                            y2_min = int (j*(self.x_size - self.overlap) + auxliar2)
   
                            cont2 = cont2 + 1
                           
                            if cont1 >= (size_r_x/2):
                                cont1 = 0
                            if cont2 >= (size_r_y/2 - 1):
                                cont2 = 0
                    else:
                        x2_min = i*(self.x_size - self.overlap)
                        y2_min = j*(self.y_size - self.overlap)
                    x1_max = x1_min + self.x_size
                    y1_max = y1_min + self.y_size
                    x2_max = x2_min + self.x_size
                    y2_max = y2_min + self.y_size
                   
                    if x1_min > self.tam_x:
                        x1_min = int (i*(self.x_size - self.overlap))
                        x1_max = x1_min + self.x_size
                    if x2_min > self.tam_x:
                        x2_min = int (i*(self.x_size - self.overlap))
                        x2_max = x2_min + self.x_size
                    if x1_max > self.tam_x:
                        x1_min = int (i*(self.x_size - self.overlap))
                        x1_max = x1_min + self.x_size
                    if x2_max > self.tam_x:
                        x2_min = int (i*(self.x_size - self.overlap))
                        x2_max = x2_min + self.x_size
       
                    if y1_min >= self.tam_y:
                        y1_min = int (j*(self.y_size - self.overlap))
                        y1_max = y1_min + self.y_size
                    if y2_min >= self.tam_y:
                        y2_min = int (j*(self.y_size - self.overlap))
                        y2_max = y2_min + self.y_size
                    if y1_max >= self.tam_y:
                        y1_min = int (j*(self.y_size - self.overlap))
                        y1_max = y1_min + self.y_size
                    if y2_max >= self.tam_y:
                        y2_min = int (j*(self.y_size - self.overlap))
                        y2_max = y2_min + self.y_size
                   
                    c = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
                    d = np.copy (self.im2[x2_min:x2_max, y2_min:y2_max])
                    
                    if c.size == 0:
                        x1_min = i* (self.x_size - self.overlap)
                        y1_min = j* (self.x_size - self.overlap)
                        x1_max = x1_min + self.x_size
                        y1_max = y1_min + self.y_size
                        c = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
                    if d.size == 0:
                        x1_min = i* (self.x_size - self.overlap)
                        y1_min = j* (self.x_size - self.overlap)
                        x1_max = x1_min + self.x_size
                        y1_max = y1_min + self.y_size
                        d = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
                        
                    janela = interrogation_window (c , d)
                   
                    janela.normxcorr2 ('same')

                    xx, yy = janela.find_peak()
                    xx, yy = janela.gauss_subpixel_peak_position (xx, yy)
                    xx -= (self.x_size + self.x_size - 1)//2
                    yy -= (self.y_size + self.y_size - 1)//2
                       
                    dpx[i, j], dpy[i, j] = -xx, yy
            dp_anterior.dpx = np.zeros((size_r_x, size_r_y))
            dp_anterior.dpy = np.zeros((size_r_x, size_r_y))
            dp_anterior.dpx = dpx
            dp_anterior.dpy = dpy
            cont1 = cont1 +1
            cont2 = 0
            aux = aux +1
        return dp_anterior
