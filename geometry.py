#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Author: Jiazhi Li
# Date: 30 Apr 2016
# Version: 1.10
# Function: 2D Monte Carlo

# Class - used to build the geometry of the model.

import numpy as np

class Geometry:

    def __init__(self, name):
        
        self.name = name;


    def circle(self, radius, pitch, pin_num):

        tot_num = pin_num ** 2;
        
        x_max = (pin_num-1) / 2;
        x_min = -x_max;
        
        cor_seq = np.linspace(x_min,x_max,pin_num) * pitch;

        
        self.x = np.tile(cor_seq,pin_num);
        self.y = np.repeat(-cor_seq,pin_num);
        self.r = np.repeat(radius,tot_num);
            
