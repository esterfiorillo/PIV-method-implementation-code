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
#     ./interface/interface.py
#
# Description
#     indicates that folder has modules for python.
#
#------------------------------------------------------------------------------

import sys
import os


"""
   try to load interface modules
"""
try:
    from PyQt5 import uic
#    from PyQt5.QtWidgets import QApplication
    from PyQt5 import QtWidgets
    from PyQt5.QtGui import QDesktopServices
    from PyQt5.QtCore import QUrl
except ModuleNotFoundError:
    print("PyQt5 missing. Please install module.")
    sys.exit()    


"""
   try to load fundamental modules
"""
try:
    import numpy as np
    import concurrent.futures
    import threading
except ModuleNotFoundError as E:
    print("Missing critical module. Please install module:")
    print(E)
    sys.exit()    


"""
    try to load pre processing functions
"""
try:
#    from pre_processing.pre_processing import calc_background 
#    from pre_processing.pre_processing import calc_m 
    from pre_processing.pre_processing import thread_mmao 
    from pre_processing.pre_processing import thread_bckg 
#    from pre_processing.image_par import image_par 
except Exception as E:
    print("Missing critical pre processing module. Please install module:")
    print(E)
    sys.exit()    


"""
    try to load processing functions
"""
try:
    from processing.processing import thread_processing 
#    from processing.displacement_map import displacement_map 
#    from processing.interrogation_window import interrogation_window
except Exception as E:
    print("Missing critical processing module. Please install module:")
    print(E)
    sys.exit()     

    
Form, Window = uic.loadUiType("./interface/interface5.4.ui")

dpx_list1 = []
dpy_list1 = []
dpx_list2 = []
dpy_list2 = []

class mywindow(QtWidgets.QMainWindow):


    def __init__(self):
        super(mywindow, self).__init__()
        self.ui = Window()
        self.ui.form = Form()
        self.ui.form.setupUi(self)
        self.ui.form.Ok_button2.clicked.connect(self.populate_im_lines)
        self.ui.form.image1_button.clicked.connect(self.select_images1)
        self.ui.form.actionUML_Diagram.triggered.connect(self.load_pdf)
       
    def select_images1(self):
        #Function that takes the file address of the first image
        
        photo_path1, ext1 = QtWidgets.QFileDialog.getOpenFileName(self, "Select Image, format xxxx0000.xxx, last numbers will define sequence")
        if photo_path1:
            self.ui.form.image1.setText(photo_path1)
        self.string_photo = photo_path1

        
    def populate_im_lines(self):
        
        #Take the value entered by the user for the number of images
        
        num_images = self.ui.form.Number_of_images.text()
        num_images = int (num_images)
        photo_path2 = self.string_photo
        num_images2 = int (num_images/2)
        
        #Get the string that has the name of the directory with the images, the prefix of the file and the number of the first image
        print(photo_path2)
        
        """ Alteration of file load """
        spl = os.path.split(photo_path2)
        dir = spl[0]
        
        e = np.log10(num_images).round()
        n = int(5 + e)
        
        file_prefix = spl[1][:-n]
        
        num_primeira = spl[1][-n]
        
        file_form = spl[1][-4:]
        
        # dir = photo_path2.split(os.path.sep)[-2].split("/")[0]
        # file_prefix = dir + "_"
        # num_primeira = photo_path2.split(os.path.sep)[-1].split("/")[0]
        # num_primeira = num_primeira.split(os.path.sep)[-1].split(".")[0]
        # num_primeira = num_primeira.split(os.path.sep)[-1].split("_")[3]
        # num_primeira = int (num_primeira)
        
        w_size = self.ui.form.WindowSize.text() #window size
        w_size = int(w_size)
        if w_size <16:
            QtWidgets.QMessageBox.about (self, "Not valid", "Type another value for Window Size")
            return
        ovl = self.ui.form.overlap.text()
        ovl = int (ovl)
        if ovl > w_size:
            QtWidgets.QMessageBox.about (self, "Not valid", "Type another value for Overlap")
            return
        n_iterations = self.ui.form.number_iterations2.text()
        n_iterations = int (n_iterations)
        met = self.ui.form.write_method.currentText()
        
        
        with concurrent.futures.ThreadPoolExecutor() as executer:
            f1 = executer.submit(thread_mmao, num_images2, dir, file_prefix, num_primeira, file_form)
            f2 = executer.submit(thread_bckg, num_images, dir, file_prefix, num_primeira, file_form)
            mao = f1.result()
            bck_ground = f2.result()  
        x = threading.Thread(target = thread_processing, args = (num_images, dir, file_prefix, num_primeira, file_form, bck_ground, mao, met, w_size, ovl, n_iterations))
        x.start()

    def load_pdf(self):

        """
        This function loads the pdf located at /docs and shows it to the user using
        the system's default pdf reader
        """
        
        QDesktopServices.openUrl(QUrl("docs/UML_diagram.pdf", mode=QUrl.TolerantMode));
    
