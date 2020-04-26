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
#     ./processing/multigrid_method.py
#------------------------------------------------------------------------------

"""
   load basic modules
"""
from processing.interrogation_window import interrogation_window
from processing.methods.methods import Methods
from processing.methods.normal_method import normal_method
import sys

"""
   try to load fundamental modules
"""
try:
    import numpy as np
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit() 


class multigrid_method (normal_method):
    """
    This class represents the second mehtod of findind the displacement vector map.
    It uses the vector map found in the normal method, so this class is a child class of normal_method class.
    The multigrid method is an iterative method, in which the size of the interrogation windows is reduced throughout the process. In addition, the interrogation windows are shifted in the images according to the displacement map found in the previous iteration.
    
    Attributes
    ----------
    im1: 2d np.array
        Image 1
    im2: 2d np.array
        Image 2
    tam_x, tam_y: int
        Images dimensions
    x_size: int
        Interrogation window size in x axis
    y_size: int
        Interrogation window size in y axis
    ovl: int
        Size of the overlap
    n_passadas: int
        Number of iterations of the method
    
    Methods
    -------
    resize_to_input_shape (map_a, dpx_map, dpy_map)
    multigrid_method1()
    
    """
    def __init__ (self, im1, im2, w_size, ovl, n_passadas):
        """
        Constructor
        
        Parameters
        ---------
        im1 : 2d np.array
            Image 1 of the par of images
        im2 : 2d np.array
             Image 2 of the par of images
        w_size: int
            Size of the interrogation window
        ovl: int
            Size of the overlap
        n_passadas: int
            
        Raises
        ------
        .
        
        """
        normal_method.__init__(self, im1, im2, w_size, ovl)
        self.num_passadas = n_passadas
        
        
    def resize_to_input_shape (self, map_a, dpx_map, dpy_map):
        """
        Function that takes an object of the type displacement_map and doubles its dimensions. It will be used to add the displacements of the previous vector map.
        
        Parameters
        ----------
        map_a : diplacement_map object
            Input displacement vector map
        dpx_map: 2d np.array
            Displacement vector map in x axis
        dpy_map: 2d np.array
            Displacement vector map in y axis
        
        Raises
        ------
        dpx_map: 2d np.array
            Output displacement vector map in x axis
        dpy_map: 2d np.array
            Output displacement vector map in y axis 
        
        """
        c = np.zeros((2, 2))
        d = np.zeros((2, 2))
        tam_ax1 = len (map_a.dpx)
        tam_ax2 = len (map_a.dpx[0])
        c_list1 = np.zeros((tam_ax1*tam_ax2, 2, 2))
        c_list2 = np.zeros((tam_ax1*tam_ax2, 2, 2))
        cont = 0
    
        for i in range (tam_ax1):  
            for j in range (tam_ax2):
                aux1 = map_a.dpx[i][j]
                aux2 = map_a.dpy[i][j]
    
                c[0][0] = aux1
                d[0][0] = aux2
                c_list1[cont] = c
                c_list2[cont] = d
                cont = cont +1
    
        cont = 0
        for i in range(0, len (dpx_map), 2):
            for j in range (0, len(dpx_map[0]), 2):
                dpx_map[i:i+2, j:j+2] = c_list1[cont]
                dpy_map[i:i+2, j:j+2] = c_list2[cont]
                cont = cont +1
    
        return dpx_map, dpy_map


    def multigrid_method1 (self):
        """
        Function that implements the Multigrid method.
        Reference:
        F. Scarano & M. L. Riethmuller, Iterative multigrid approach in PIV image processing with discrete window offset, Experiments in Fluids 26 (1999)
        
        Parameters
        ----------
        .
        
        Raises
        ------
        dp_anterior: displacement_map object
            displacement vector map resulting from the method
        
        """
    
        aux = 0
        dp_anterior = self.first_iteration()
        while (aux != self.num_passadas):
    
            self.x_size = self.x_size//2
            self.y_size = self.y_size//2
            self.overlap = self.overlap//2
            size_r_x, size_r_y = self.result_dimensions(self.x_size, self.y_size, self.overlap)
            
            dpx = np.zeros((size_r_x, size_r_y))
            dpy = np.zeros((size_r_x, size_r_y))
            c = np.zeros((self.x_size, self.y_size))
            d = np.zeros((self.x_size, self.y_size))
            
            tam_dpxa = 2*len(dp_anterior.dpx)
            tam_dpya = 2*len(dp_anterior.dpx[0])
            
            dp_anterior.resize_dp ()
            dpx_aux = np.zeros((size_r_x, size_r_y))
            dpy_aux = np.zeros((size_r_x, size_r_y))
            
            dpx_add = np.zeros((tam_dpxa, tam_dpya))
            dpy_add = np.zeros((tam_dpxa, tam_dpya))
            dpx_add, dpy_soma = self.resize_to_input_shape(dp_anterior, dpx_add, dpy_add)
            
 
            #dpx_aux will always have bigger dimensions than dpx_add, so this loop is not a problem
            for i in range (tam_dpxa):
                for j in range (tam_dpya):
                    if (np.isnan(dpx_add[i][j])):
                       dpx_aux[i][j] = 0
                    else:
                        dpx_aux[i][j] = dpx_add[i][j]
                    if (np.isnan(dpy_soma[i][j])):
                       dpy_aux[i][j] = 0
                    else:
                        dpy_aux[i][j] = dpy_add[i][j]
    
    
            for i in range (0, size_r_x):
                for j in range (0, size_r_y):
                    x1_min = i* (self.x_size - self.overlap)
                    y1_min = j* (self.x_size - self.overlap)

                    x2_min = int (i* (self.x_size - self.overlap) + dpx_aux[i][j])
                    y2_min = int (j* (self.x_size - self.overlap) + dpy_aux[i][j])
                    
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
            aux = aux +1
    
        return dp_anterior

#Second implementaton:
        
#    def multigrid_method1 (self):
#    #Function that implements the Multigrid method, which is an iterative procedure for traversing images with interrogation windows while applying a cross correlation.
#    #It consists of moving windows according to the displacement map found in the previous iteration.
#    #At the same time, the size of the interrogation windows is reduced between one iteration and another.
#    #In addition, a data validation criterion is applied to the displacement maps.
#    #In:
#    #    num_passadas: int
#    #        number of iterations
#    #    x_size: int
#    #        interrogation window size in x axis
#    #    y_size: int
#    #        interrogation window size in y axis
#    #    overlap: int
#    #    tam_x: int
#    #        x-dimension of the images
#    #    tam_y: int
#    #        y-dimension of the images
#    #    dpx_anterior: 2d array
#    #        previous displacement in x axis found in the function fisrt_iteration
#    #    dpy_anterior: 2d array
#    #        previous displacement in y axis found in the function fisrt_iteration
#    #    im1: 2d array
#    #        first image
#    #    im2: 2d array
#    #        second image
#    #Out:
#    #    dpx_anterior: 2d array
#    #        displacement in x
#    #    dpy_anterior: 2d array
#    #        displacement in y
#    
#    #Reference:
#    #F. Scarano & M. L. Riethmuller, Iterative multigrid approach in PIV image processing with discrete window offset, Experiments in Fluids 26 (1999)
#        aux = 0
#        while (aux != self.num_passadas):
#            dp_anterior = self.first_iteration()
#            self.x_size = self.x_size//2
#            self.y_size = self.y_size//2
#            self.overlap = self.overlap//2
#            size_r_x, size_r_y = self.result_dimensions()
#            ind1_x, ind1_y, ind2_x, ind2_y = dp_anterior.index_guard()
#   
#            dp_anterior.resize_dp ()
#
#            dpx = np.zeros((size_r_x, size_r_y))
#            dpy = np.zeros((size_r_x, size_r_y))
#            c = np.zeros((size_r_x, size_r_y))
#            d = np.zeros((size_r_x, size_r_y))
#            cont1 = 0
#            cont2 = 0
#            for i in range (0, size_r_x):
#                for j in range (0, size_r_y):
#                    x1_min = i* (self.x_size - self.overlap)
#                    y1_min = j* (self.x_size - self.overlap)
#                    if i == ind2_x[cont1]:
#                        if j == ind2_y[cont2]:
#                            teste1_indice1 = int (ind1_x[cont1])
#                            teste1_indice2 = int (ind1_y[cont2])
#   
#                            if (np.isnan(dp_anterior.dpx[teste1_indice1][teste1_indice2])):
#                                auxliar1 = 0
#                            else:
#                                auxliar1 = dp_anterior.dpx[teste1_indice1][teste1_indice2]
#                            if (np.isnan(dp_anterior.dpy[teste1_indice1][teste1_indice2])):
#                                auxliar2 = 0
#                            else:
#                                auxliar2 = dp_anterior.dpy[teste1_indice1][teste1_indice2]
#                            x2_min = int (i*(self.x_size - self.overlap) + auxliar1)
#                            y2_min = int (j*(self.x_size - self.overlap) + auxliar2)
#   
#                            cont2 = cont2 + 1
#                           
#                            if cont1 >= (size_r_x/2):
#                                cont1 = 0
#                            if cont2 >= (size_r_y/2 - 1):
#                                cont2 = 0
#                    else:
#                        x2_min = i*(self.x_size - self.overlap)
#                        y2_min = j*(self.y_size - self.overlap)
#                    x1_max = x1_min + self.x_size
#                    y1_max = y1_min + self.y_size
#                    x2_max = x2_min + self.x_size
#                    y2_max = y2_min + self.y_size
#                   
#                    if x1_min > self.tam_x:
#                        x1_min = int (i*(self.x_size - self.overlap))
#                        x1_max = x1_min + self.x_size
#                    if x2_min > self.tam_x:
#                        x2_min = int (i*(self.x_size - self.overlap))
#                        x2_max = x2_min + self.x_size
#                    if x1_max > self.tam_x:
#                        x1_min = int (i*(self.x_size - self.overlap))
#                        x1_max = x1_min + self.x_size
#                    if x2_max > self.tam_x:
#                        x2_min = int (i*(self.x_size - self.overlap))
#                        x2_max = x2_min + self.x_size
#       
#                    if y1_min >= self.tam_y:
#                        y1_min = int (j*(self.y_size - self.overlap))
#                        y1_max = y1_min + self.y_size
#                    if y2_min >= self.tam_y:
#                        y2_min = int (j*(self.y_size - self.overlap))
#                        y2_max = y2_min + self.y_size
#                    if y1_max >= self.tam_y:
#                        y1_min = int (j*(self.y_size - self.overlap))
#                        y1_max = y1_min + self.y_size
#                    if y2_max >= self.tam_y:
#                        y2_min = int (j*(self.y_size - self.overlap))
#                        y2_max = y2_min + self.y_size
#                   
#                    c = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
#                    d = np.copy (self.im2[x2_min:x2_max, y2_min:y2_max])
#                    
#                    if c.size == 0:
#                        x1_min = i* (self.x_size - self.overlap)
#                        y1_min = j* (self.x_size - self.overlap)
#                        x1_max = x1_min + self.x_size
#                        y1_max = y1_min + self.y_size
#                        c = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
#                    if d.size == 0:
#                        x1_min = i* (self.x_size - self.overlap)
#                        y1_min = j* (self.x_size - self.overlap)
#                        x1_max = x1_min + self.x_size
#                        y1_max = y1_min + self.y_size
#                        d = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
#                        
#                    janela = interrogation_window (c , d)
#                   
#                    janela.normxcorr2 ('same')
#
#                    xx, yy = janela.find_peak()
#                    xx, yy = janela.gauss_subpixel_peak_position (xx, yy)
#                    xx -= (self.x_size + self.x_size - 1)//2
#                    yy -= (self.y_size + self.y_size - 1)//2
#                       
#                    dpx[i, j], dpy[i, j] = -xx, yy
#            dp_anterior.dpx = np.zeros((size_r_x, size_r_y))
#            dp_anterior.dpy = np.zeros((size_r_x, size_r_y))
#            dp_anterior.dpx = dpx
#            dp_anterior.dpy = dpy
#            cont1 = cont1 +1
#            cont2 = 0
#            aux = aux +1
#        return dp_anterior
