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
#         clahe_histogram_equalization  
#         remove_background
#         calc_teta
#         homogenize_brightness
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
    import cv2
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit()    
    

class pre_processing:
    """  
    This is a class to represent the pre processing methods applied to the images
    
        
       Attributes
       --------
       ...
       
       Methods
       -------
       sobel_filter1
       sobel_filter2
       laplacian_filter
       clahe_histogram_equalization  
       remove_background
       calc_teta
       homogenize_brightness 
    """    
    
    def __init__(self, mmao, bck_g):
        """
        Constructor
        
        Parameters
        ---------
        mmao : float
            Necessary value for the brightness homogenization
        bck_g : float
            Value for the background of the images
        
        Raises
        ------
        .
        
        """
        self.mmao = mmao
        self.bckg = bck_g
 
    
    def sobel_filter1(self, im):
        """
        Function that applies a Sobel filter to an image. This filter has the function
        of detecting horizontal edges.
        
        Parameters
        ----------
        im: 2d np.array
            Input image
        
        Raises
        ------
        im: 2d np.array
            Input image with a sobel filter
        
        """
        im = cv2.Sobel(im, cv2.CV_64F,1,0,ksize=5)

        return im
  
    
    def sobel_filter2(self, im):
        """
        Function that applies a Sobel filter to an image. This filter has the function
        of detecting vertical edges.
        
        Parameters
        ----------
        im: 2d np.array
            Input image
        
        Raises
        ------
        im: 2d array
            Input image with a sobel filter
        
        """
        im = cv2.Sobel(im, cv2.CV_64F,0,1,ksize=5)
        return im
        
    def laplacian_filter(self, im):
        """
        Function that applies a Laplacian filter to an image. This filter has the function
        of detecting horizontal and vertical edges.
        
        Parameters
        ----------
        im: 2d np.array
            Input image
        
        Raises
        ------
        im: 2d np.array
            Input image with a laplacian filter
        
        """
        im = cv2.Laplacian(im,cv2.CV_64F)
        return im
    
    def clahe_histogram_equalization(self, im):
        """
        Function that applies a CLAHE (contrast limited adaptative histogram equalization) to an image. 
        
        Parameters
        ----------
        im: 2d np.array
            Input image
        
        Raises
        ------
        im: 2d np.array
            Input image after applied CLAHE 
        
        """
        clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
        im = clahe.apply(im)
        return im
 
        
    def remove_background (self, im):
        """
        Function that subtracts the background from the image
        
        Parameters
        ----------
        im: 2d np.array
            Input image
        
        Raises
        ------
        im: 2d np.array
            Input image after having its background extracted
        """
        im = im - self.bckg
        return im
        
    def calc_teta(self, im2):
        """
        This function and homogenize_brightness function aim to homogenize the brightness of the second image according to that of the first.
        It was made to be apliyed to the second image of the par, not both of then.
        Reference : Brightness correction, homogenization and adjunstment of images for face recognition, Eduardo Machado Silva, UNESP, ISSN 2316-9664, Volume 14, fev. 2019, Edic ̧ao Ermac
        
        Parameters
        ----------
        im2: 2d np.array
            Input image
        
        Raises
        ------
        teta: float
            Value that will be used in homogenize_brightness function
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
        
        Parameters
        ----------
        im2: 2d np.array
            Input image
        
        Raises
        ------
        im2: 2d np.array
            Input image after brightness correction
        """
        teta = self.calc_teta(im2)
        tam_x = len(im2)
        tam_y = len(im2[0])
        for i in range(tam_x):
            for j in range (tam_y):
                im2[i][j] = im2[i][j] - teta*im2[i][j]
        return im2


 