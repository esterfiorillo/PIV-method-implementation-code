#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 10:50:39 2020

@author: esterfiorillo
"""

import numpy as np
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
from numpy import log

#If the user does not have the opencv library installed,
#the modo_erro_cv variable will be equal to 1, and otherwise, it will be equal to 0. 
#This will be used to pass a message on the graphical interface
#if the library is not installed
modo_erro_cv = 0
try:
    import cv2
except ModuleNotFoundError:
    modo_erro_cv = 1

def calc_background (n_im, dir, file_prefix, num_primeira):
#Function that calculates the average of all images to be processed
    im1 = plt.imread(dir + "/" + file_prefix + str (num_primeira).zfill(5) + ".tif")
    im1 = np.asarray (im1)
    if len (np.shape (im1)) > 2:
        im1 = np.mean(im1, -1) #rgb to gray  
    soma = im1
    for i in range (0, n_im):            
        num_primeira = num_primeira + 1
        im2 = plt.imread(dir + "/" + file_prefix + str (num_primeira).zfill(5) + ".tif")
        im2 = np.asarray (im2)

        if len (np.shape (im2)) > 2:
            im2 = np.mean (im2, -1) #rgb to gray
        soma = soma + im2
    res = soma/n_im
    return res

def calc_m (n_im, dir, file_prefix, num_primeira):
    aux = 0
    for i in range (0, n_im):
        im1 = plt.imread(dir + "/" + file_prefix + str (num_primeira).zfill(5) + ".tif")
        im1 = np.asarray (im1)
        if len (np.shape (im1)) > 2:
            im1 = np.mean(im1, -1) #rgb to gray
        h,bins = np.histogram(im1.ravel(),256,[0,256])
        h_media = np.mean(h)
        aux = aux + h_media
        num_primeira = num_primeira + 2
    mzao = aux/n_im
    return mzao

class image_par:
#Base class of normal_method class and multigrid_method class,
#which represents the methods of traversing the pair of images
#with the interrogation windows while the cross correlation between them is applyed.
#The advantage of using inheritance in this case is the possibility
#of adding more methods without having to change anything in the code
#Atributes:
#    image1: 2d np.narray
#    image2: 2d np.narray
#Functions:
#    sobel_filter1
#    sobel_filter2
#    laplacian_filter
#    remove_background
#    calc_teta
#    homogenize_brightness
    
    def __init__(self, image1, image2): #Constructor
        self.im1 = image1
        self.im2 = image2
        self.tam_x, self.tam_y = np.shape(self.im1)
   
    def sobel_filter1(self):
#Function that applies a Sobel filter to an image. This filter has the function
#of detecting horizontal edges.
        self.im1 = cv2.Sobel(self.im1, cv2.CV_64F,1,0,ksize=5)
        self.im2 = cv2.Sobel(self.im2, cv2.CV_64F,1,0,ksize=5)
        #return self.im1
   
    def sobel_filter2(self):
#Function that applies a Sobel filter to an image. This filter has the function
#of detecting vertical edges.
        self.im1 = cv2.Sobel(self.im1, cv2.CV_64F,0,1,ksize=5)
        self.im2 = cv2.Sobel(self.im2, cv2.CV_64F,0,1,ksize=5)
       
    def laplacian_filter(self):
#Function that applies a Laplacian filter to an image. This filter has the function
#of detecting vertical and horizontal edges.
        self.im1 = cv2.Laplacian(self.im1,cv2.CV_64F)
        self.im2 = cv2.Laplacian(self.im2,cv2.CV_64F)
       
    def remove_background (self, res):
#Function that subtracts that average from the pair of images.
#This is intended to remove the background from the images.
        self.im1 = self.im1 - res
        self.im2 = self.im2 - res
        
#The two functions below (calc_teta and homogenize_brightness) aim to 
#homogenize the brightness of the second image according to that of the first
#Brightness correction, homogenization and adjunstment of images for face recognition, Eduardo Machado Silva, UNESP, ISSN 2316-9664, Volume 14, fev. 2019, Edic ̧ao Ermac  
   
    def calc_teta(self, mmao):
        h,bins = np.histogram(self.im2.ravel(),256,[0,256])
        h_min = min(h)
        h_max = max(h)
        h_media = np.mean(h)
        teta = (mmao - h_media)/(h_max - h_min)
        return teta
   
    def homogenize_brightness(self,  mmao):
        teta = self.calc_teta(mmao)
        tam_x = len(self.im2)
        tam_y = len(self.im2[0])
        for i in range(tam_x):
            for j in range (tam_y):
                self.im2[i][j] = self.im2[i][j] - teta*self.im2[i][j]
       
class interrogation_window:
#Class that represents the interrogation windows and the functions applyed to then.
#Atributes:
    #window1: 2d np_array
    #window2: 2d np_array
    #cross_corr: 2d np_array
#Functions:
    #normxcorr2()
    #find_peak()
    #gauss_subpixel_peak_position
    def __init__(self, w1, w2): #Constructor
        self.window1 = w1
        self.window2 = w2
   
    def normxcorr2(self, mode="full"):
########################################################################################
# Author: Ujash Joshi, University of Toronto, 2017                                     #
# Based on Octave implementation by: Benjamin Eltzner, 2014 <b.eltzner@gmx.de>         #
# Octave/Matlab normxcorr2 implementation in python 3.5                                #
# Details:                                                                             #
# Normalized cross-correlation. Similiar results upto 3 significant digits.            #
# https://github.com/Sabrewarrior/normxcorr2-python/master/norxcorr2.py                #
# http://lordsabre.blogspot.ca/2017/09/matlab-normxcorr2-implemented-in-python.html    #
########################################################################################
   

#    Input arrays should be floating point numbers.
#    :param template: N-D array, of template or filter you are using for cross-correlation.
#    Must be less or equal dimensions to image.
#    Length of each dimension must be less than length of image.
#    :param image: N-D array
#    :param mode: Options, "full", "valid", "same"
#    full (Default): The output of fftconvolve is the full discrete linear convolution of the inputs.
#    Output size will be image size + 1/2 template size in each dimension.
#    valid: The output consists only of those elements that do not rely on the zero-padding.
#    same: The output is the same size as image, centered with respect to the ‘full’ output.
#    :return: N-D array of same dimensions as image. Size depends on mode parameter.
   
    # If this happens, it is probably a mistake
        if np.ndim(self.window1) > np.ndim(self.window2) or \
                len([i for i in range(np.ndim(self.window1)) if self.window1.shape[i] > self.window2.shape[i]]) > 0:
                print("normxcorr2: self.window1 larger than IMG. Arguments may be swapped.")

        self.window1 = self.window1 - np.mean(self.window1)
        self.window2 = self.window2 - np.mean(self.window2)

        a1 = np.ones(self.window1.shape)
        # Faster to flip up down and left right then use fftconvolve instead of scipy's correlate
        ar = np.flipud(np.fliplr(self.window1))
        out = fftconvolve(self.window2, ar.conj(), mode=mode)
   
        self.window2 = fftconvolve(np.square(self.window2), a1, mode=mode) - \
                np.square(fftconvolve(self.window2, a1, mode=mode)) / (np.prod(self.window1.shape))

    # Remove small machine precision errors after subtraction
        self.window2[np.where(self.window2 < 0)] = 0

        self.window1 = np.sum(np.square(self.window1))
        out = out / np.sqrt(self.window2 * self.window1)

        # Remove any divisions by 0 or very close to 0
        out[np.where(np.logical_not(np.isfinite(out)))] = 0
   
        self.cross_corr = out
       
    def find_peak (self):    
#Function that returns the coordinates in x and y of the biggest point (peak)
#in the 2d array resulted from the cross correlation.
        ind = self.cross_corr.argmax()
        s = self.cross_corr.shape[1]
        xx_peak = ind // s
        yy_peak = ind % s
        return xx_peak, yy_peak

    def gauss_subpixel_peak_position (self, x_pico, y_pico):
#Function that returns the gauss interpolation coordinates
#
#In:
#    x_pico: int
#         x-coordinate of the peak
#    y_pico: int
#         y-coordinate of the peak
#Out:
#    subp_peak_position[0]: int
#        x-coordinate in the cross-correlation matrix resulting from gauss interpolation
#    subp_peak_position[1]: int
#        y-coordinate in the cross-correlation matrix resulting from gauss interpolation
        if x_pico + 1 >= len(self.cross_corr):
            subp_peak_position = x_pico, y_pico
        elif y_pico +1 >= len(self.cross_corr[0]):
            subp_peak_position = x_pico, y_pico
        else:
           
            c = self.cross_corr[x_pico, y_pico]
            cl = self.cross_corr[x_pico - 1, y_pico]
            cr = self.cross_corr[x_pico + 1, y_pico]
            cd = self.cross_corr[x_pico, y_pico - 1]
            cu = self.cross_corr[x_pico, y_pico + 1]
            if np.any(np.array([c, cl, cr, cd, cu]) < 0):
                subp_peak_position = (((x_pico - 1) * cl + x_pico * c + (x_pico + 1) * cr) / (cl + c + cr),
((y_pico - 1) * cd + y_pico * c + (y_pico + 1) * cu) / (cd + c + cu))
            else:    
                subp_peak_position = (x_pico + ((log(cl) - log(cr)) / (2 * log(cl) - 4 * log(c) + 2 * log(cr))),
y_pico + ((log(cd) - log(cu)) / (2 * log(cd) - 4 * log(c) + 2 * log(cu))))
        return subp_peak_position[0], subp_peak_position[1]        
   

class displacement_map:
#Class that represents the displacement map resulting from the methods.
#It also have the functions applied to the map for the data validation.
#Atributes:
    #dp_x: 2d np.narray
    #dp_y: 2d np.narray
#Functions:
    #neighborhood_median(int ii, int jj)
    #calculate_rm (int ii, int jj)
    #calculate_r (int ii, int jj)
    #neighborhood_median(int ii, int jj)
    #replacement3()
    #resize_dp
    #index_guard()
   
    def __init__ (self, dp_x, dp_y): #Constructor
        self.dpx = dp_x
        self.dpy = dp_y
       
    def neighborhood_median (self, ii, jj):
    #Function that calculates the median of the neighboorhood points (3x3) of a certain position (ii, jj) on the displacement map
        aux1 = np.array([[self.dpx[ii-1][jj-1], self.dpx[ii-1][jj], self.dpx[ii-1][jj+1], self.dpx[ii][jj-1], self.dpx[ii][jj+1], self.dpx[ii+1][jj-1], self.dpx[ii+1][jj], self.dpx[ii+1][jj+1]]])
        aux2 = np.array([[self.dpy[ii-1][jj-1], self.dpy[ii-1][jj], self.dpy[ii-1][jj+1], self.dpy[ii][jj-1], self.dpy[ii][jj+1], self.dpy[ii+1][jj-1], self.dpy[ii+1][jj], self.dpy[ii+1][jj+1]]])
        med1 = np.median(aux1)
        med2 = np.median(aux2)
        return med1, med2

    def calculate_rm (self, ii, jj):
    #Function that calculates the median of the subtration of the neighborhood points (3x3) and the neighboorhood median
        median1, median2 = self.neighborhood_median(ii, jj)
        aux1 = np.array([self.dpx[ii-1][jj-1] - median1, self.dpx[ii-1][jj]- median1, self.dpx[ii-1][jj+1]- median1, self.dpx[ii][jj-1]- median1, self.dpx[ii][jj+1]- median1, self.dpx[ii+1][jj-1]- median1, self.dpx[ii+1][jj]- median1, self.dpx[ii+1][jj+1]- median1])
        aux2 = np.array([self.dpy[ii-1][jj-1] - median2, self.dpy[ii-1][jj]- median2, self.dpy[ii-1][jj+1]- median2, self.dpy[ii][jj-1]- median2, self.dpy[ii][jj+1]- median2, self.dpy[ii+1][jj-1]- median2, self.dpy[ii+1][jj]- median2, self.dpy[ii+1][jj+1]- median2])
        rm1 = np.median (aux1)
        rm2 = np.median (aux2)
        return rm1, rm2

    def calculate_r (self, ii, jj):
    #Function that calculates: (point of the displacement map - median of the 3x3 neighborhood) / rm - 0.1
        median1, median2 = self.neighborhood_median(ii, jj)
        rm1, rm2 = self.calculate_rm (ii, jj)
        r1 = (self.dpx[ii][jj] - median1)/(rm1 - 0.1)
        r2 = (self.dpy[ii][jj] - median2)/(rm2 - 0.1)
        return r1, r2

    def neighborhood_mean (self, ii, jj):
    #Function that calculates the mean of the neighboorhood points (3x3) of a certain position (ii, jj) on the displacement map
        aux1 = np.array([self.dpx[ii-1][jj-1], self.dpx[ii-1][jj], self.dpx[ii-1][jj+1], self.dpx[ii][jj-1], self.dpx[ii][jj+1], self.dpx[ii+1][jj-1], self.dpx[ii+1][jj], self.dpx[ii+1][jj+1]])
        aux2 = np.array([self.dpy[ii-1][jj-1], self.dpy[ii-1][jj], self.dpy[ii-1][jj+1], self.dpy[ii][jj-1], self.dpy[ii][jj+1], self.dpy[ii+1][jj-1], self.dpy[ii+1][jj], self.dpy[ii+1][jj+1]])
        med1 = np.mean(aux1)
        med2 = np.mean(aux2)
        return med1, med2

    def replacement3(self):
    #Function that calculates the value of r for all points on the displacement map
    #and replaces the value of points with the neighborhood mean, if the value of r
    #at this point is less than 2.
#Reference:
#F. Scarano & J. Westerweel, Universal Outlier detection for PIV data , Experiments in Fluids 39 (2005)
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
#Function that divides the values ​​of the displacement vectors by two
        self.dpx = self.dpx/2
        self.dpy = self.dpy/2

    def index_guard (self):
#function that stores in a vector twice the value of the indexes of 
#the displacement map of the previous iteration
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

class normal_method(image_par):
#Method1
    def __init__ (self, im1, im2, w_size, ovl): #Constructor
        image_par.__init__(self, im1, im2)
        self.x_size = w_size
        self.y_size = w_size
        self.overlap = ovl

    def result_dimensions (self):  
#Function used to find the dimensions of the result of applying cross correlation in the interrogation windows, according to the interrogation window size and overlap
#In:
#    x_im: int
#        x-dimension of the image
#    y_im: int
#        y-dimension of the image
#    window_size: int
#        size of the interrogation window that will be used to go through the images
#    overlap: int
#        
#Out:
#    x_dimension: int
#        x-dimension of matrix resulting from applying cross-correlation between interrogation windows
#    y_dimension: int
#       y-dimension of matrix resulting from applying cross-correlation between interrogation windows
        x_dimension = (self.tam_x - self.x_size) // (self.x_size - self.overlap) + 1
        y_dimension = (self.tam_y - self.y_size) // (self.y_size - self.overlap) + 1
        return x_dimension, y_dimension

    def first_iteration (self):
#Function in which the pair of images is traversed with the interrogation windows while applying the cross correlation between them.
#In this case, the windows will be shifted normally, with no method to make the process more efficient.
#This step is required to find the first displacement map that will be used in Multigrid.
#In:
#    im1: 2d array
#        first image
#    im2: 2d array
#        second image
#    tam_x: int
#        x-dimension of the images
#    tam_y: int
#        y-dimension of the images
#    x_size: int
#        interrogation window size in x axis
#    y_size: int
#        interrogation window size in y axis
#    overlap: int
#Out:
#    dpx1: 2d array
#        displacement in x
#    dpy1: 2d array
#        displacement in y    
        size_r_x, size_r_y = self.result_dimensions()
        #create matrix for substitute divisions of im1
        a = np.zeros ((self.x_size, self.y_size))
        #create matrix for substitute search window in im2
        b = np.zeros ((self.x_size, self.y_size))

        dpx1 = np.zeros ((size_r_x, size_r_y))
        dpy1 = np.zeros ((size_r_x, size_r_y))
       
        for i in range (0, size_r_x):
            for j in range (0, size_r_y):
                x1_min = i*(self.x_size - self.overlap)
                x1_max = x1_min + self.x_size

                y1_min = j*(self.y_size - self.overlap)
                y1_max = y1_min + self.y_size

                a = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
                b = np.copy (self.im2[x1_min:x1_max, y1_min:y1_max])
               
                janela = interrogation_window (a , b)
                janela.normxcorr2 ('same')
                
                xx, yy = janela.find_peak()
                xx, yy = janela.gauss_subpixel_peak_position (xx, yy)
                xx -= (self.x_size + self.x_size - 1)//2
                yy -= (self.y_size + self.y_size - 1)//2

                dpx1[i][j], dpy1[i][j] = -xx, yy

        res = displacement_map(dpx1, dpy1)
        return res

class multigrid_method (image_par):
#Method2
    def __init__ (self, im1, im2, w_size, ovl, n_passadas): #Constructor
        image_par.__init__(self, im1, im2)
        self.x_size = w_size
        self.y_size = w_size
        self.overlap = ovl
        self.num_passadas = n_passadas

    def result_dimensions (self):
#Function used to find the dimensions of the result of applying cross correlation in the interrogation windows, according to the interrogation window size and overlap.
#In:
#    x_im: int
#        x-dimension of the image
#    y_im: int
#        y-dimension of the image
#    window_size: int
#        size of the interrogation window that will be used to go through the images
#    overlap: int
#        
#Out:
#    x_dimension: int
#        x-dimension of matrix resulting from applying cross-correlation between interrogation windows
#    y_dimension: int
#       y-dimension of matrix resulting from applying cross-correlation between interrogation windows
   
        x_dimension = (self.tam_x - self.x_size) // (self.x_size - self.overlap) + 1
        y_dimension = (self.tam_y - self.y_size) // (self.y_size - self.overlap) + 1
        return x_dimension, y_dimension
   
    def multigrid_method1 (self, dp_anterior):
#Function that implements the Multigrid method, which is an iterative procedure for traversing images with interrogation windows while applying a cross correlation.
#It consists of moving windows according to the displacement map found in the previous iteration.
#At the same time, the size of the interrogation windows is reduced between one iteration and another.
#In addition, a data validation criterion is applied to the displacement maps.
#In:
#    num_passadas: int
#        number of iterations
#    x_size: int
#        interrogation window size in x axis
#    y_size: int
#        interrogation window size in y axis
#    overlap: int
#    tam_x: int
#        x-dimension of the images
#    tam_y: int
#        y-dimension of the images
#    dpx_anterior: 2d array
#        previous displacement in x axis found in the function fisrt_iteration
#    dpy_anterior: 2d array
#        previous displacement in y axis found in the function fisrt_iteration
#    im1: 2d array
#        first image
#    im2: 2d array
#        second image
#Out:
#    dpx_anterior: 2d array
#        displacement in x
#    dpy_anterior: 2d array
#        displacement in y

#Reference:
#F. Scarano & M. L. Riethmuller, Iterative multigrid approach in PIV image processing with discrete window offset, Experiments in Fluids 26 (1999)
        aux = 0
        while (aux != self.num_passadas):
            self.x_size = self.x_size//2
            self.y_size = self.y_size//2
            self.overlap = self.overlap//2
            size_r_x, size_r_y = self.result_dimensions()
            ind1_x, ind1_y, ind2_x, ind2_y = dp_anterior.index_guard()
   
            dp_anterior.resize_dp ()

            dpx = np.zeros((size_r_x, size_r_y))
            dpy = np.zeros((size_r_x, size_r_y))
            c = np.zeros((size_r_x, size_r_y))
            d = np.zeros((size_r_x, size_r_y))
            cont1 = 0
            cont2 = 0
            for i in range (0, size_r_x):
                for j in range (0, size_r_y):
                    x1_min = i* (self.x_size - self.overlap)
                    y1_min = j* (self.x_size - self.overlap)
                    if i == ind2_x[cont1]:
                        if j == ind2_y[cont2]:
                            teste1_indice1 = int (ind1_x[cont1])
                            teste1_indice2 = int (ind1_y[cont2])
   
                            if (np.isnan(dp_anterior.dpx[teste1_indice1][teste1_indice2])):
                                auxliar1 = 0
                            else:
                                auxliar1 = dp_anterior.dpx[teste1_indice1][teste1_indice2]
                            if (np.isnan(dp_anterior.dpy[teste1_indice1][teste1_indice2])):
                                auxliar2 = 0
                            else:
                                auxliar2 = dp_anterior.dpy[teste1_indice1][teste1_indice2]
                            x2_min = int (i*(self.x_size - self.overlap) + auxliar1)
                            y2_min = int (j*(self.x_size - self.overlap) + auxliar2)
   
                            cont2 = cont2 + 1
                           
                            if cont1 >= (size_r_x/2):
                                cont1 = 0
                            if cont2 >= (size_r_y/2 - 1):
                                cont2 = 0
                    else:
                        x2_min = i*(self.x_size - self.overlap)
                        y2_min = j*(self.y_size - self.overlap)
                    x1_max = x1_min + self.x_size
                    y1_max = y1_min + self.y_size
                    x2_max = x2_min + self.x_size
                    y2_max = y2_min + self.y_size
                   
                    if x1_min > self.tam_x:
                        x1_min = int (i*(self.x_size - self.overlap))
                        x1_max = x1_min + self.x_size
                    if x2_min > self.tam_x:
                        x2_min = int (i*(self.x_size - self.overlap))
                        x2_max = x2_min + self.x_size
                    if x1_max > self.tam_x:
                        x1_min = int (i*(self.x_size - self.overlap))
                        x1_max = x1_min + self.x_size
                    if x2_max > self.tam_x:
                        x2_min = int (i*(self.x_size - self.overlap))
                        x2_max = x2_min + self.x_size
       
                    if y1_min >= self.tam_y:
                        y1_min = int (j*(self.y_size - self.overlap))
                        y1_max = y1_min + self.y_size
                    if y2_min >= self.tam_y:
                        y2_min = int (j*(self.y_size - self.overlap))
                        y2_max = y2_min + self.y_size
                    if y1_max >= self.tam_y:
                        y1_min = int (j*(self.y_size - self.overlap))
                        y1_max = y1_min + self.y_size
                    if y2_max >= self.tam_y:
                        y2_min = int (j*(self.y_size - self.overlap))
                        y2_max = y2_min + self.y_size
                   
                    c = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
                    d = np.copy (self.im2[x2_min:x2_max, y2_min:y2_max])
                    
                    if c.size == 0:
                        x1_min = i* (self.x_size - self.overlap)
                        y1_min = j* (self.x_size - self.overlap)
                        x1_max = x1_min + self.x_size
                        y1_max = y1_min + self.y_size
                        c = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
                    if d.size == 0:
                        x1_min = i* (self.x_size - self.overlap)
                        y1_min = j* (self.x_size - self.overlap)
                        x1_max = x1_min + self.x_size
                        y1_max = y1_min + self.y_size
                        d = np.copy (self.im1[x1_min:x1_max, y1_min:y1_max])
                        
                    janela = interrogation_window (c , d)
                   
                    janela.normxcorr2 ('same')

                    xx, yy = janela.find_peak()
                    xx, yy = janela.gauss_subpixel_peak_position (xx, yy)
                    xx -= (self.x_size + self.x_size - 1)//2
                    yy -= (self.y_size + self.y_size - 1)//2
                       
                    dpx[i, j], dpy[i, j] = -xx, yy
            dp_anterior.dpx = np.zeros((size_r_x, size_r_y))
            dp_anterior.dpy = np.zeros((size_r_x, size_r_y))
            dp_anterior.dpx = dpx
            dp_anterior.dpy = dpy
            cont1 = cont1 +1
            cont2 = 0
            aux = aux +1
        return dp_anterior
