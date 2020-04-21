#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 11:37:19 2020

@author: esterfiorillo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

entrada1 = pd.read_csv('dpx_resultado_multigrid', sep = '\t')
a = entrada1.as_matrix()
a = np.array(a)

entrada2 = pd.read_csv('dpy_resultado_multigrid', sep = '\t')
b = entrada2.as_matrix()
b = np.array(b)

plt.quiver(a, b)