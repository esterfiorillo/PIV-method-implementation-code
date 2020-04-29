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
#     ./processing/image_par.py
#-----------------------------------------------------------------------------

"""
This file represents the class that has the pre processing functions, which aims to make corrections and improve the quality of the images
"""

"""
   load basic modules
"""
import sys

"""
   try to load fundamental modules
"""

try:
    import numpy as np
    import cv2
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit()    
    

class pre_processing:
    """  
    This is a class to represent the pre processing methods applied to the images

    """    
    
    def __init__(self, mmao, bck_g):
        """
        Constructor
        
        :type mmao : float
        :param mmao: Necessary value for the brightness homogenization
        
        :type bck_g : float
        :param bck_g: Value for the background of the images
        
        """
        self.mmao = mmao
        self.bckg = bck_g
 
    
    def sobel_filter1(self, im):
        """
        Function that applies a Sobel filter to an image. This filter has the function
        of detecting horizontal edges.
        
        :type im: 2d np.array
        :param im: Input image
        
        :rtype: 2d np.array
        :return: Input image with a sobel filter
        
        """
        im = cv2.Sobel(im, cv2.CV_64F,1,0,ksize=5)

        return im
  
    
    def sobel_filter2(self, im):
        """
        Function that applies a Sobel filter to an image. This filter has the function
        of detecting vertical edges.
        
        :type im: 2d np.array
        :param im: Input image
        
        :rtype: 2d np.array
        :return: Input image with a sobel filter
        
        """
        im = cv2.Sobel(im, cv2.CV_64F,0,1,ksize=5)
        return im
        
    def laplacian_filter(self, im):
        """
        Function that applies a Laplacian filter to an image. This filter has the function
        of detecting horizontal and vertical edges.
        
        :type im: 2d np.array
        :param im: Input image
        
        :rtype: 2d np.array
        :return: Input image with a laplacian filter
        
        """
        im = cv2.Laplacian(im,cv2.CV_64F)
        return im
    
    def clahe_histogram_equalization(self, im):
        """
        Function that applies a CLAHE (contrast limited adaptative histogram equalization) to an image. 
        
        :type im: 2d np.array
        :param im: Input image
        
        :rtype: 2d np.array
        :return: Input image after CLAHE been apliyed
        
        """
        clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
        im = clahe.apply(im)
        return im
 
        
    def remove_background (self, im):
        """
        Function that subtracts the background from the image
        
        :type im: 2d np.array
        :param im: Input image
        
        :rtype: 2d np.array
        :return: Input image after having its background extracted
        """
        im = im - self.bckg
        return im
        
    def calc_teta(self, im2):
        """
        This function and homogenize_brightness function aim to homogenize the brightness of the second image according to that of the first.
        It was made to be apliyed to the second image of the par, not both of then.
        Reference : Brightness correction, homogenization and adjunstment of images for face recognition, Eduardo Machado Silva, UNESP, ISSN 2316-9664, Volume 14, fev. 2019, Edic ̧ao Ermac
        
        :type im: 2d np.array
        :param im: Input image
        
        :rtype: float
        :return: Value that will be used in homogenize_brightness function
        """
        h,bins = np.histogram(im2.ravel(),256,[0,256])
        h_min = min(h)
        h_max = max(h)
        h_media = np.mean(h)
        teta = (self.mmao - h_media)/(h_max - h_min)
        return teta
 
    
    def homogenize_brightness(self, im2):
        """
        This function and homogenize_brightness function aim to homogenize the brightness of the second image according to that of the first.
        It was made to be apliyed to the second image of the par, not both of then.
        Reference : Brightness correction, homogenization and adjunstment of images for face recognition, Eduardo Machado Silva, UNESP, ISSN 2316-9664, Volume 14, fev. 2019, Edic ̧ao Ermac
        
        :type im: 2d np.array
        :param im: Input image
        
        :rtype: 2d np.array
        :return: Input image after brightness correction
        """
        teta = self.calc_teta(im2)
        tam_x = len(im2)
        tam_y = len(im2[0])
        for i in range(tam_x):
            for j in range (tam_y):
                im2[i][j] = im2[i][j] - teta*im2[i][j]
        return im2


 