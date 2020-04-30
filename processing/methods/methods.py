# -*- coding: utf-8 -*-
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
#     ./processing/Methods/methods.py

"""
This file contains the base class for the processing methods.
"""


"""
   load basic modules
"""

from pre_processing.pre_processing_class import pre_processing
from processing.interrogation_window import interrogation_window
import sys
from processing.displacement_map import displacement_map

"""
   try to load fundamental modules
"""
try:
    import numpy as np
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit() 
    
class Methods():
    """
    Base class of normal_method class and multigrid_method class, which represents the methods of traversing the pair of images with the interrogation windows while the cross correlation between them is applyed.   
    Since this is the base class of the methods, it stores the pair of images being processed.
    
    """
    
    def __init__(self, image1, image2): 
        """
        Constructor
        
        :type image1 : 2d np.array
        :param image1: Image 1 of the par of images
        
        :type image2 : 2d np.array
        :param image2: Image 2 of the par of images
            
        """
        
        self.im1 = image1
        self.im2 = image2
        self.tam_x, self.tam_y = np.shape(self.im1)
    
        
    def result_dimensions (self, x_size, y_size, overlap):
        """
        Function used to find the dimensions of the result of applying cross correlation in the interrogation windows, according to the interrogation window size and overlap.
        
        :type x_size: int
        param x_size: Size of the interrogation window that will be used to go through the images in x axis
        :type y_size: int 
        :param y_size: Size of the interrogation window that will be used to go through the images in y axis
        :type overlap: int
        :param overlap: How much overlap
        
        :rtype: int
        :return: Dimensions of matrix resulting from applying cross-correlation between interrogation windows
        """
        x_dimension = (self.tam_x - x_size) // (x_size - overlap) + 1
        y_dimension = (self.tam_y - y_size) // (y_size - overlap) + 1
        return x_dimension, y_dimension
        
    
        
        