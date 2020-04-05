#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 10:26:17 2020

@author: esterfiorillo
"""

from teste_quarentena2 import calc_background
from teste_quarentena2 import calc_m
from teste_quarentena2 import image_par
from teste_quarentena2 import normal_method
from teste_quarentena2 import multigrid_method 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

#If the user does not have the opencv library installed,
#the modo_erro_cv variable will be equal to 1, and otherwise, it will be equal to 0. 
#This will be used to pass a message on the graphical interface
#if the library is not installed
modo_erro_cv = 0
try:
    import cv2
except ModuleNotFoundError:
    modo_erro_cv = 1
    
#librarys from interface
from PyQt5 import uic
from PyQt5.QtWidgets import QApplication
from PyQt5 import QtWidgets
import sys


#Load interface
Form, Window = uic.loadUiType("interface5.4.ui")

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
       
    def select_images1(self):
        #Function that takes the file address of the first image
        photo_path1, ext1 = QtWidgets.QFileDialog.getOpenFileName(self, "Select Photo")
        if photo_path1:
            self.ui.form.image1.setText(photo_path1)
        self.string_photo = photo_path1
       
    def populate_im_lines(self):
        if modo_erro_cv == 1:
            QtWidgets.QMessageBox.about (self, "Not valid", "Install OpenCV library")
        
        #Take the value entered by the user for the number of images
        num_images = self.ui.form.Number_of_images.text()
        num_images = int (num_images)
        photo_path2 = self.string_photo
        num_images2 = int (num_images/2)
        
        #Get the string that has the name of the directory with the images, the prefix of the file and the number of the first image
        dir = photo_path2.split(os.path.sep)[-2].split("/")[0]
        file_prefix = dir + "_"
        num_primeira = photo_path2.split(os.path.sep)[-1].split("/")[0]
        num_primeira = num_primeira.split(os.path.sep)[-1].split(".")[0]
        num_primeira = num_primeira.split(os.path.sep)[-1].split("_")[3]
        num_primeira = int (num_primeira)
        
    
        bck_ground = calc_background(num_images, dir, file_prefix, num_primeira)
        mao = calc_m(num_images2, dir, file_prefix, num_primeira)
        
        for i in range (1, num_images):
            
            im1 = plt.imread(dir + "/" + file_prefix + str (num_primeira).zfill(5) + ".tif")
            print (file_prefix + str (num_primeira).zfill(5) + ".tif")
            im1 = np.asarray (im1)
            num_primeira = num_primeira + 1
            im2 = plt.imread(dir + "/" + file_prefix + str (num_primeira).zfill(5) + ".tif")
            print (file_prefix + str (num_primeira).zfill(5) + ".tif")
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
   
            w_size = self.ui.form.WindowSize.text() #window size
            w_size = int(w_size)
            if w_size <16:
                QtWidgets.QMessageBox.about (self, "Not valid", "Type another value for Window Size")
                return
            if w_size > tam_x:
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
                    pd.DataFrame(dpx_csv).to_csv("dpx_resultado_multigrid", sep='\t')
                    dpy_csv = pd.DataFrame(dpy_def)
                    pd.DataFrame(dpy_csv).to_csv("dpy_resultado_multigrid", sep='\t')
                    print ("Process Finished")
               
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

    
app = QApplication([])
application = mywindow()
application.show()
sys.exit(app.exec())