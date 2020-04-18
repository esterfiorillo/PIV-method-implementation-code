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
#     ./pre_processing/pre_processing.py
#
# Description
#     Function that pre processes the images for PIV.
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
    import matplotlib.pyplot as plt
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit()    
    
    


def thread_mmao(num_images2, dir, file_prefix, num_primeira, file_form):
    mao = calc_m(num_images2, dir, file_prefix, num_primeira, file_form)
    print("finish calculating mmao")
    return mao


def thread_bckg(num_images, dir, file_prefix, num_primeira, file_form):
    bck_ground = calc_background(num_images, dir, file_prefix, num_primeira, file_form)
    print ("finish calculating bckg")
    return bck_ground


def calc_background (n_im, dir, file_prefix, num_primeira, file_form):
    
    # Order of magnitude of the number of images
    e = int(np.log10(n_im).round())
    
    #Function that calculates the average of all images to be processed
    
    for i in range (0, n_im):            
        num_primeira = int(num_primeira) + i
        im2 = plt.imread(dir + "/" + file_prefix + str (num_primeira).zfill(e) + file_form)
        im2 = np.asarray (im2)

        if len (np.shape (im2)) > 2:
            im2 = np.mean (im2, -1) #rgb to gray
        
        if i == 0:
            soma = im2
        else:
            soma = soma + im2
    res = soma/n_im
    return res


def calc_m (n_im, dir, file_prefix, num_primeira, file_form):
    
    aux = 0
    
    # Order of magnitude of the number of images
    e = int(np.log10(n_im).round())
    
    for i in range (0, n_im):
        num_primeira = int(num_primeira) + i
        im1 = plt.imread(dir + "/" + file_prefix + str (num_primeira).zfill(e) + file_form)
        im1 = np.asarray (im1)
        if len (np.shape (im1)) > 2:
            im1 = np.mean(im1, -1) #rgb to gray
        h,bins = np.histogram(im1.ravel(),256,[0,256])
        h_media = np.mean(h)
        aux = aux + h_media
        
    mzao = aux/n_im
    return mzao
