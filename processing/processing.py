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
#     ./processing/processing.py
#
# Description
#     Function that processes the image pair.
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
    import pandas as pd
    import time

except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit()    
 

"""
    try to load pre processing functions
"""
try:
    from pre_processing.image_par import image_par 
except Exception as E:
    print("Missing critical pre processing module. Please install module:")
    print(E)
    sys.exit()   


    
"""
    try to load processing methods
"""

from processing.methods.normal_method import normal_method
from processing.methods.multigrid_method import multigrid_method

# Function:

   

def thread_processing (num_images, dir, file_prefix, num_primeira, file_form, bck_ground, mao, met, w_size, ovl, n_iterations):
    start = time.perf_counter()
    dpx_list1 = []
    dpy_list1 = []
    dpx_list2 = []
    dpy_list2 = [] 
    # Order of magnitude of the number of images
    e = int(np.log10(num_images).round())
    
    for i in range (1, num_images):
        im1 = plt.imread(dir + "/" + file_prefix + str (num_primeira).zfill(e) + file_form)
        print (f'processing {file_prefix + str (num_primeira).zfill(5) + file_form}')
        im1 = np.asarray (im1)
        num_primeira = int(num_primeira) + 1
        im2 = plt.imread(dir + "/" + file_prefix + str (num_primeira).zfill(e) + file_form)
        print (f'processing {file_prefix + str (num_primeira).zfill(5) + file_form}')
        im2 = np.asarray (im2)
        num_primeira = num_primeira + 1
       
        if len (np.shape (im1)) > 2:
            im1 = np.mean(im1, -1) #rgb to gray
        if len (np.shape (im2)) > 2:
            im2 = np.mean (im2, -1) #rgb to gray
       
        pr = image_par(im1, im2)
       
        pr.sobel_filter1 ()
        pr.sobel_filter1 ()

        pr.remove_background(bck_ground)
        pr.homogenize_brightness(mao)
        tam_x, tam_y = np.shape(im1)
   
        if met == 'Multigrid':
            par1 = normal_method(pr.im1, pr.im2, w_size, ovl)
            r = par1.first_iteration()
           
            par2 = multigrid_method(im1, im2, w_size, ovl, n_iterations)
            s = par2.multigrid_method1(r)
            s.replacement3()
            dpx_list1.append(s.dpx)
            dpy_list1.append(s.dpy)
            if i == (num_images -1):
                dpx_def = sum (dpx_list1)/len (dpx_list1)
                dpy_def = sum (dpy_list1)/len (dpx_list1)
                 
                #saves the average of dpx and dpy results in two csv files
                dpx_csv = pd.DataFrame(dpx_def)
                pd.DataFrame(dpx_csv).to_csv("dpx_resultado_multigrid", sep='\t', header = False, index = False)
                dpy_csv = pd.DataFrame(dpy_def)
                pd.DataFrame(dpy_csv).to_csv("dpy_resultado_multigrid", sep='\t', header = False, index = False)
                print ("Process Finished")
                finish = time.perf_counter()
                sec = round (finish - start, 2)
                print(f"execution time = {sec} seconds") 
           
        if met == 'Normal':
            par1 = normal_method(pr.im1, pr.im2, w_size, ovl)
            r = par1.first_iteration()
            dpx_list2.append(s.dpx)
            dpy_list2.append(s.dpy)
            
            if i == (num_images -1):
                dpx_def2 = sum (dpx_list2)/len (dpx_list2)
                dpy_def2 = sum (dpy_list2)/len (dpx_list2)
            
                #saves the average of dpx and dpy results in two csv files
                dpx_csv2 = pd.DataFrame(dpx_def2)
                pd.DataFrame(dpx_csv2).to_csv("dpx_resultado_multigrid", sep='\t')
                dpy_csv2 = pd.DataFrame(dpy_def2)
                pd.DataFrame(dpy_csv2).to_csv("dpy_resultado_multigrid", sep='\t')
                print ("Process Finished")
                finish = time.perf_counter()
                sec = round (finish - start, 2)
                print(f"execution time = {sec} seconds") 
