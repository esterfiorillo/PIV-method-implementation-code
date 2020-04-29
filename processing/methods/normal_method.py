## -*- coding: utf-8 -*-
#"""
#
# (               )     )   |
# )\ )  *   )  ( /(  ( /(   |
#(()/(` )  /(  )\()) )\())  |
# /(_))( )(_))((_)\ ((_)\   | LTHN PIV: Is an opensource PIV toolbox
#(_)) (_(_())  _((_) _((_)  |
#| |  |_   _| | || || \| |  |
#| |__  | |   | __ || .` |  | Website: https://github.com/esterfiorillo/PIV-method-implementation-code
#|(___|(|_|   |_||_||_|\_|  |
#                           | 
# )\ ) )\ )                 | CDTN - Centro de Desenvolvimento da Tecnologia Nuclear
#(()/((()/( (   (           | LTHN - Laboratório de Termo-Hidráulica e Neutrônica
# /(_))/(_)))\  )\          | Belo Horizonte, MG, Brasil
#(_)) (_)) ((_)((_)         |
#| _ \|_ _|\ \ / /          | @authors: esterfiorillo, acampagnole 
#|  _/ | |  \ V /           |
#|_|  |___|  \_/            |
#
#
#"""

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
"""
This file contains the class that represents the implementation of the normal method
"""

"""
   load basic modules
"""
from processing.interrogation_window import interrogation_window
import sys
from processing.displacement_map import displacement_map

from processing.methods.methods import Methods

"""
   try to load fundamental modules
"""
try:
    import numpy as np
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit() 




class normal_method(Methods):
    """
    Class that represents the first method of traversing the pair of images with the interrogation windows while the cross correlation between them is applyed.
    In this case, the interrogation windows are moved with a constant displacement both horizontally and vertically, and they do not vary in size throughout the process.
    This is a child class of the methods class.
    
    """
    
    def __init__ (self, im1, im2, w_size, ovl):
        """
        Constructor
        
        :type im1 : 2d np.array
        :param im1: Image 1 of the par of images
        :type im2 : 2d np.array
        :param im2: Image 2 of the par of images
        :type w_size: int
        :param w_size: Size of the interrogation window
        :type ovl: int
        :param ovl: Size of the overlap
        """
        Methods.__init__(self, im1, im2)
        self.x_size = w_size
        self.y_size = w_size
        self.overlap = ovl


    def first_iteration (self):
        """
        Function that represents the normal method.
        The normal method can be seen as a required step to find the first displacement vector map that will be used in Multigrid method.
        
        :rtype: displacement_map
        :return: Displacement vector map resulting from the method
      
        """
        
        size_r_x, size_r_y = self.result_dimensions(self.x_size, self.y_size, self.overlap)
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
