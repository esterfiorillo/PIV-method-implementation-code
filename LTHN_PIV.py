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
#     ./LTHN_PIV.py
#
# Description
#     Core application of LTHN PIV. 
#     Imports modules and classes from folders.
#     Initializes the interface.
#  
#------------------------------------------------------------------------------
"""
This file contains the function that opens the interface. The user have to execute only it
"""


"""
   load basic modules
"""
import sys



import numpy as np
from numpy import log
import matplotlib.pyplot as plt
import pandas as pd
import time
import concurrent.futures
import threading
from scipy.signal import fftconvolve
import cv2

from PyQt5.QtWidgets import QApplication   
from PyQt5 import uic
from PyQt5.QtWidgets import QApplication
from PyQt5 import QtWidgets
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtCore import QUrl


from pre_processing.pre_processing import calc_background 
from pre_processing.pre_processing import calc_m 
from pre_processing.pre_processing import thread_mmao 
from pre_processing.pre_processing import thread_bckg 
from pre_processing.pre_processing_class import pre_processing

from processing.processing import thread_processing 
from processing.displacement_map import displacement_map 
from processing.interrogation_window import interrogation_window

from processing.methods.normal_method import normal_method
from processing.methods.multigrid_method import multigrid_method  



from interface.interface import mywindow as mywindow



def main():
    """
    
    Starts LTHN PIV interface
    
    """
    app = QApplication([])
    application = mywindow()
    application.show()
    app.exec()

    
if __name__ == "__main__":
    main()
    
