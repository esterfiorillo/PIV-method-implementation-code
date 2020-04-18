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
#     ./processing/displacement_map.py
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

"""
   load basic modules
"""
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
    

class displacement_map:

   
    def __init__ (self, dp_x, dp_y): #Constructor
        self.dpx = dp_x
        self.dpy = dp_y
     
        
    def neighborhood_median (self, ii, jj):
    # Function that calculates the median of the neighboorhood points (3x3) of a certain position (ii, jj) on the displacement map
    
        aux1 = np.array([[self.dpx[ii-1][jj-1], self.dpx[ii-1][jj], self.dpx[ii-1][jj+1], self.dpx[ii][jj-1], self.dpx[ii][jj+1], self.dpx[ii+1][jj-1], self.dpx[ii+1][jj], self.dpx[ii+1][jj+1]]])
        aux2 = np.array([[self.dpy[ii-1][jj-1], self.dpy[ii-1][jj], self.dpy[ii-1][jj+1], self.dpy[ii][jj-1], self.dpy[ii][jj+1], self.dpy[ii+1][jj-1], self.dpy[ii+1][jj], self.dpy[ii+1][jj+1]]])
        med1 = np.median(aux1)
        med2 = np.median(aux2)
        return med1, med2


    def calculate_rm (self, ii, jj):
    # Function that calculates the median of the subtration of the neighborhood points (3x3) and the neighboorhood median
    
        median1, median2 = self.neighborhood_median(ii, jj)
        aux1 = np.array([self.dpx[ii-1][jj-1] - median1, self.dpx[ii-1][jj]- median1, self.dpx[ii-1][jj+1]- median1, self.dpx[ii][jj-1]- median1, self.dpx[ii][jj+1]- median1, self.dpx[ii+1][jj-1]- median1, self.dpx[ii+1][jj]- median1, self.dpx[ii+1][jj+1]- median1])
        aux2 = np.array([self.dpy[ii-1][jj-1] - median2, self.dpy[ii-1][jj]- median2, self.dpy[ii-1][jj+1]- median2, self.dpy[ii][jj-1]- median2, self.dpy[ii][jj+1]- median2, self.dpy[ii+1][jj-1]- median2, self.dpy[ii+1][jj]- median2, self.dpy[ii+1][jj+1]- median2])
        rm1 = np.median (aux1)
        rm2 = np.median (aux2)
        return rm1, rm2


    def calculate_r (self, ii, jj):
    # Function that calculates: (point of the displacement map - median of the 3x3 neighborhood) / rm - 0.1
    
        median1, median2 = self.neighborhood_median(ii, jj)
        rm1, rm2 = self.calculate_rm (ii, jj)
        r1 = (self.dpx[ii][jj] - median1)/(rm1 - 0.1)
        r2 = (self.dpy[ii][jj] - median2)/(rm2 - 0.1)
        return r1, r2


    def neighborhood_mean (self, ii, jj):
    # Function that calculates the mean of the neighboorhood points (3x3) of a certain position (ii, jj) on the displacement map
    
        aux1 = np.array([self.dpx[ii-1][jj-1], self.dpx[ii-1][jj], self.dpx[ii-1][jj+1], self.dpx[ii][jj-1], self.dpx[ii][jj+1], self.dpx[ii+1][jj-1], self.dpx[ii+1][jj], self.dpx[ii+1][jj+1]])
        aux2 = np.array([self.dpy[ii-1][jj-1], self.dpy[ii-1][jj], self.dpy[ii-1][jj+1], self.dpy[ii][jj-1], self.dpy[ii][jj+1], self.dpy[ii+1][jj-1], self.dpy[ii+1][jj], self.dpy[ii+1][jj+1]])
        med1 = np.mean(aux1)
        med2 = np.mean(aux2)
        return med1, med2


    def replacement3(self):
    # Function that calculates the value of r for all points on the displacement map
    # and replaces the value of points with the neighborhood mean, if the value of r
    # at this point is less than 2.
    # Reference:
    # F. Scarano & J. Westerweel, Universal Outlier detection for PIV data , Experiments in Fluids 39 (2005)
    
        ux_size = len (self.dpx)
        uy_size = len (self.dpx[0])
        for i in range (1, ux_size-1):
            for j in range (1, uy_size-1):
                r1_value, r2_value = self.calculate_r(i, j)
                if r1_value <=2:
                    self.dpx[i][j], al = self.neighborhood_mean (i, j)

        ux_size = len (self.dpy)
        uy_size = len (self.dpy[0])
        for i in range (1, ux_size-1):
            for j in range (1, uy_size-1):
                r1_value, r2_value = self.calculate_r(i, j)
                if r2_value <=2:
                    al, self.dpy[i][j] = self.neighborhood_mean (i, j)


    def resize_dp (self):
    # Function that divides the values ​​of the displacement vectors by two
    
        self.dpx = self.dpx/2
        self.dpy = self.dpy/2


    def index_guard (self):
    # function that stores in a vector twice the value of the indexes of 
    # the displacement map of the previous iteration
    
        indices1_x = np.zeros((len(self.dpx)))
        indices1_y = np.zeros((len(self.dpx[0])))
        indices2_x = np.zeros((len(self.dpx)))
        indices2_y = np.zeros((len(self.dpx[0])))
        indices1_x = list (range(0, len(self.dpx)))
        indices1_y = list (range(0, len(self.dpx[0])))
        for i in range (len(self.dpx)):
            indices2_x[i] = 2*indices1_x[i]
        for i in range (len(self.dpx[0])):
            indices2_y[i] = 2*indices1_y[i]
        return indices1_x, indices1_y, indices2_x, indices2_y
