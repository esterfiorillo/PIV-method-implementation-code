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
#     ./commum/iofunc.py
#
# Description
#     indicates that folder has modules for python.
#
#------------------------------------------------------------------------------

"""
   load basic modules
"""
import sys
import os

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

class io_func:
    """Class to deal with the IO such as read_image, write_image, read_input. Expects an working directory path on instancing"""
    
    def __init__ (self, **kwargs): #Constructor
        self.path = kwargs.get('path')
        print(self.path)
        if self.path is None:
            self.nopath_error() 


    def nopath_error(self):
        print('No path defined. Exiting.')
        

    def read_image(self, **kwargs):
        """Read an image function. Expects path, file_prefix, nzeros, imnumber, file_form."""
        self.file_prefix = kwargs.get('file_prefix')
        self.nzeros = kwargs.get('nzeros')
        self.imnumber = kwargs.get('imnumber')
        self.file_form = kwargs.get('file_form')
        
        print(os.path.join(self.path,str(self.file_prefix) + str(self.imnumber).zfill(int(self.nzeros)) + str(self.file_form)))

        return plt.imread(os.path.join(self.path,str(self.file_prefix) + str(self.imnumber).zfill(int(self.nzeros)) + str(self.file_form)))


photo_path2 = os.path.abspath("image_test1.png")
spl = os.path.split(photo_path2)
path = spl[0]

num_images = 2

nzeros = np.log10(num_images).round()
n = int(5 + nzeros)

file_prefix = spl[1][:-n]

imnumber = spl[1][-n]

file_form = spl[1][-4:]

read = io_func(path = path)
im = read.read_image( file_prefix = file_prefix, imnumber = imnumber, file_form = file_form, nzeros = nzeros)

view = plt.imshow(im)
plt.show()
