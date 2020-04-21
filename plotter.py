#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 11:37:19 2020

@author: esterfiorillo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ws = input("Type the final interrogation window size:")
ws = int (ws)

entrada1 = pd.read_csv('dpx_resultado_multigrid', sep = '\t')
a = entrada1.as_matrix()
a = np.array(a)

entrada2 = pd.read_csv('dpy_resultado_multigrid', sep = '\t')
b = entrada2.as_matrix()
b = np.array(b)

tam_y, tam_x= np.shape(a)
x = np.arange(0, tam_x, 1)
y = np.arange(0, tam_y, 1)
x, y= np.meshgrid(x, y)

def reescaling(x, y, a, b, ws):
    x = x*ws
    y= y*ws
    a = a*ws
    b = b*ws
    return a, b, x, y

a, b, x, y = reescaling(x, y, a, b, ws)

plt.quiver(x, y, a, b)