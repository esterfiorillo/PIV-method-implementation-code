#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 17:17:48 2020

@author: esterfiorillo
"""

#from pre_processing.pre_processing_class import pre_processing
#from processing.interrogation_window import interrogation_window
import sys
#from processing.displacement_map import displacement_map

try:
    import numpy as np
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit() 
    
class Methods():
    
    def __init__(self, image1, image2): #Constructor
        self.im1 = image1
        self.im2 = image2
        self.tam_x, self.tam_y = np.shape(self.im1)
    
        #self.final_map = 
        
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
        
    
        
        