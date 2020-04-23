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
#     ./LTHN_PIV.py
#
# Description
#     Core application of LTHN PIV. 
#     Imports modules and classes from folders.
#     Initializes the interface.
#  
#------------------------------------------------------------------------------

"""
   load basic modules
"""
import sys
# import os


"""
    try to load fundamental modules
"""
try:
    import numpy as np
    from numpy import log
    import matplotlib.pyplot as plt
    import pandas as pd
    import time
    import concurrent.futures
    import threading
    from scipy.signal import fftconvolve
    import cv2
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit()    
    
    
"""
   try to load interface modules
"""
try:
    # from PyQt5 import uic
    from PyQt5.QtWidgets import QApplication

    # from PyQt5 import QtWidgets
except ModuleNotFoundError:
    print("PyQt5 missing. Please install module.")
    sys.exit()     


# """
#     try to load pre processing functions
# """
# try:
#     from pre_processing.pre_processing import calc_background 
#     from pre_processing.pre_processing import calc_m 
#     from pre_processing.pre_processing import thread_mmao 
#     from pre_processing.pre_processing import thread_bckg 
#     from pre_processing.image_par import image_par 
# except Exception as E:
#     print("Missing critical pre processing module. Please install module:")
#     print(E)
#     sys.exit()   


# """
#     try to load processing functions
# """
# try:
#     from processing.processing import thread_processing 
#     from processing.displacement_map import displacement_map 
#     from processing.interrogation_window import interrogation_window
# except Exception as E:
#     print("Missing critical processing module. Please install module:")
#     print(E)
#     sys.exit() 
    
    
# """
#     try to load processing methods
# """

# from processing.methods.normal_method import normal_method
# from processing.methods.multigrid_method import multigrid_method

# try:
#     from processing.methods.normal_method import normal_method
#     from processing.methods.multigrid_method import multigrid_method
# except Exception as E:
#     print("Missing critical processing method module. Please install module:")
#     print(E)
#     sys.exit()    


# """
#     try to load post processing functions
# """
#
# Placeholder
#
    

"""
   try to load interface
"""

# Form, Window = uic.loadUiType("./interface/interface5.4.ui")
from interface.interface import mywindow as mywindow


"""
    Start LTHN PIV interface
"""

app = QApplication([])
application = mywindow()
application.show()
sys.exit(app.exec())
