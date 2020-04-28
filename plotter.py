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
#     ./plotter.py

"""
   load basic modules
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def reescaling(x, y, a, b, ws):
    """
    Function that reescale the values x, y, a, b for the ws size
    
    Parameters
    ----------
    x : 2d np.array
        x coordinades of the vector map
    y : 2d np.array
        y coordinades of the vector map
    a: 2d np.array
        displacement vector map in x axis
    b: 2d np.array
        Displacement vector map in y axis
    ws: float
        The final size of the interrogation window after processing the images
    
    Raises
    ------
    x : 2d np.array
        x coordinades of the vector map after being reescaled
    y : 2d np.array
        y coordinades of the vector map after being reescaled
    a: 2d np.array
        displacement vector map in x axis after being reescaled
    b: 2d np.array
        Displacement vector map in y axis after being reescaled
        
    """
    x = x*ws
    y= y*ws
    a = a*ws
    b = b*ws
    return a, b, x, y

"""
Asks the user the final size of the interrogation window after processing the images
"""
ws = input("Type the final interrogation window size:")
ws = int (ws)

"""
Opens the two csv files, which have the displacement vector maps found in the processing and transforms them into np.arrays
"""
entrada1 = pd.read_csv('dpx_resultado_multigrid', sep = '\t')
a = entrada1.as_matrix()
a = np.array(a)
entrada2 = pd.read_csv('dpy_resultado_multigrid', sep = '\t')
b = entrada2.as_matrix()
b = np.array(b)

"""
Finds the values of x and y, which will be the coordinates of the displacement vector maps
"""
tam_y, tam_x= np.shape(a)
x = np.arange(0, tam_x, 1)
y = np.arange(0, tam_y, 1)
x, y= np.meshgrid(x, y)

"""
The reescaling function is applied to the displacement vector maps contained in the csv files
"""
a, b, x, y = reescaling(x, y, a, b, ws)

"""
The displacement vector maps, after being in the right scale, are plotted
"""
plt.quiver(x, y, a, b)