#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Author: Jiazhi Li
# Date: 15 Apr 2016
# Version: 1.6
# Function: 2D Monte Carlo

# Class - used to build isotopes information

import numpy as np

class Isotope:

    def __init__(self, name, xs, nd):
        
        self.name = name

        self.no = xs[0];

        self.alpha = ((xs[0] - 1)/(xs[0] + 1)) ** 2;
        
        self.xs = xs[1:len(xs)]; # unit: barns
        
        self.nd = nd;
    
        self.macro_xs = np.sum(self.xs) * 1e-24 * self.nd;

        self.scat_ratio = self.xs[2] * 1e-24 * self.nd / self.macro_xs;
        self.capt_ratio = self.xs[1] * 1e-24 * self.nd / self.macro_xs;
        self.fiss_ratio = self.xs[0] * 1e-24 * self.nd / self.macro_xs;
        self.absp_ratio = self.fiss_ratio + self.capt_ratio;

