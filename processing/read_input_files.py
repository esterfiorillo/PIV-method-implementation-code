## -*- coding: utf-8 -*-
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
#     ./interface/interface.py

"""
This file is an alternative to take the input parameters, other than through the graphical interface. It does this by reading the parameters in a txt file.
"""

"""
   try to load fundamental modules
"""
try:
    import concurrent.futures
    import threading
    import numpy as np
    import os
    import sys
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit() 

"""
    try to load pre processing functions
"""
try: 
    from pre_processing.pre_processing import thread_mmao 
    from pre_processing.pre_processing import thread_bckg 

except Exception as E:
    print("Missing critical pre processing module. Please install module:")
    print(E)
    sys.exit()    


"""
    try to load processing functions
"""
try:
    from processing.processing import thread_processing 

except Exception as E:
    print("Missing critical processing module. Please install module:")
    print(E)
    sys.exit()

def processing_input_file(file_address):
    """
    Function that reads a txt file that has the input parameters desired by the user. These parameters are image address1, number of images, size of the initial interrogation window, overlap, method used and number of iterations if the multigrid method is chosen. After that, the function calls the processing threads to perform the process of finding the csv file with the resulting displacement vector map.
    
    :type file_address: str
    :param file_address: String with the address of the input file that has the parameters
    """
    
    with open(file_address) as f:
            for line in f:
                if(line.strip() != '\n'): 
                    line_spl = line.split('=') 
                    cols = line_spl[0:2]
                    if cols[0] == 'Image1':
                        im1 = str(cols[1])
                        im1 = im1[:-1]
                    if cols[0] == 'Number of Images':
                        num_images = int (cols[1])
                    if cols[0] == 'Interrogation Window Size':
                        w_size = int (cols[1])
                    if cols[0] == 'Overlap':
                        ovl = int (cols[1])
                    if cols[0] == 'Method':
                        met = str(cols[1])
                        met = met[:-1]
                    if cols[0] == 'Number of Iterations':
                        n_iterations = int(cols[1])
    
    spl = os.path.split(im1)
    dir = spl[0]
    
    e = np.log10(num_images).round()
    n = int(5 + e)
    
    file_prefix = spl[1][:-n]
    
    num_primeira = spl[1][-n]
    
    file_form = spl[1][-4:]
    
    num_images2 = int (num_images/2)
                        
    
    with concurrent.futures.ThreadPoolExecutor() as executer:
        f1 = executer.submit(thread_mmao, num_images2, dir, file_prefix, num_primeira, file_form)
        f2 = executer.submit(thread_bckg, num_images, dir, file_prefix, num_primeira, file_form)
        mao = f1.result()
        bck_ground = f2.result()  
    x = threading.Thread(target = thread_processing, args = (num_images, dir, file_prefix, num_primeira, file_form, bck_ground, mao, met, w_size, ovl, n_iterations))
    x.start()

                
                
                
                
                
                
                
                
                