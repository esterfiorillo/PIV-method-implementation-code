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
#     ./processing/image_par.py
#
# Description
#     Base class of normal_method class and multigrid_method class,
#     which represents the methods of traversing the pair of images
#     with the interrogation windows while the cross correlation between them is applyed.
#     The advantage of using inheritance in this case is the possibility
#     of adding more methods without having to change anything in the code
#     Atributes:
#         image1: 2d np.narray
#         image2: 2d np.narray
#     Functions:
#         sobel_filter1
#         sobel_filter2
#         laplacian_filter
#         remove_background
#         calc_teta
#         homogenize_brightness
#
#------------------------------------------------------------------------------
import sys


try:
    import numpy as np
    import cv2
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit()    
    

class image_par:

    
    def __init__(self, image1, image2): #Constructor
        self.im1 = image1
        self.im2 = image2
        self.tam_x, self.tam_y = np.shape(self.im1)
 
    
    def sobel_filter1(self):
    #Function that applies a Sobel filter to an image. This filter has the function
    #of detecting horizontal edges.
        self.im1 = cv2.Sobel(self.im1, cv2.CV_64F,1,0,ksize=5)
        self.im2 = cv2.Sobel(self.im2, cv2.CV_64F,1,0,ksize=5)
        #return self.im1
  
    
    def sobel_filter2(self):
    #Function that applies a Sobel filter to an image. This filter has the function
    #of detecting vertical edges.
        self.im1 = cv2.Sobel(self.im1, cv2.CV_64F,0,1,ksize=5)
        self.im2 = cv2.Sobel(self.im2, cv2.CV_64F,0,1,ksize=5)

        
    def laplacian_filter(self):
    #Function that applies a Laplacian filter to an image. This filter has the function
    #of detecting vertical and horizontal edges.
        self.im1 = cv2.Laplacian(self.im1,cv2.CV_64F)
        self.im2 = cv2.Laplacian(self.im2,cv2.CV_64F)
 
        
    def remove_background (self, res):
    #Function that subtracts that average from the pair of images.
    #This is intended to remove the background from the images.
        self.im1 = self.im1 - res
        self.im2 = self.im2 - res
        
    #The two functions below (calc_teta and homogenize_brightness) aim to 
    #homogenize the brightness of the second image according to that of the first
    #Brightness correction, homogenization and adjunstment of images for face recognition, Eduardo Machado Silva, UNESP, ISSN 2316-9664, Volume 14, fev. 2019, Edic ̧ao Ermac  

   
    def calc_teta(self, mmao):
        h,bins = np.histogram(self.im2.ravel(),256,[0,256])
        h_min = min(h)
        h_max = max(h)
        h_media = np.mean(h)
        teta = (mmao - h_media)/(h_max - h_min)
        return teta
 
    
    def homogenize_brightness(self,  mmao):
        teta = self.calc_teta(mmao)
        tam_x = len(self.im2)
        tam_y = len(self.im2[0])
        for i in range(tam_x):
            for j in range (tam_y):
                self.im2[i][j] = self.im2[i][j] - teta*self.im2[i][j]
 