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
#     ./plotter.py

"""
This file is used to plot the final result on a graphic
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def reescaling(x, y, a, b, ws):
    """
    Function that reescale the values x, y, a, b for the ws size
    
    :type x: 2d np.array
    :param x: coordinades of the vector map
    
    :type y: 2d np.array
    :param y: coordinades of the vector map
    
    :type a: 2d np.array
    :param a: displacement vector map in x axis
    
    :type b: 2d np.array
    :param b: Displacement vector map in y axis
    
    :type ws: float
    :param ws: The final size of the interrogation window after processing the images
    
    :return: x, y, a, b after being reescaled
        
    """
    x = x*ws
    y= y*ws
    a = a*ws
    b = b*ws
    return a, b, x, y



def main():
    """
    Asks the user the final size of the interrogation window after processing the images, apply reescaling and plot the displacement vector map
    """
    ws = input("Type the final interrogation window size:")
    ws = int (ws)
    
    #Opens the two csv files, which have the displacement vector maps found in the processing and transforms them into np.arrays
    entrada1 = pd.read_csv('dpx_resultado_multigrid', sep = '\t')
    a = entrada1.as_matrix()
    a = np.array(a)
    entrada2 = pd.read_csv('dpy_resultado_multigrid', sep = '\t')
    b = entrada2.as_matrix()
    b = np.array(b)
    
    #Finds the values of x and y, which will be the coordinates of the displacement vector maps
    tam_y, tam_x= np.shape(a)
    x = np.arange(0, tam_x, 1)
    y = np.arange(0, tam_y, 1)
    x, y= np.meshgrid(x, y)
    
    #The reescaling function is applied to the displacement vector maps contained in the csv files
    a, b, x, y = reescaling(x, y, a, b, ws)
    
    
    #The displacement vector maps, after being in the right scale, are plotted
    plt.quiver(x, y, a, b)


if __name__ == "__main__":
    main()
