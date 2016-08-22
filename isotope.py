#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Author: Jiazhi Li
# Date: 15 Apr 2016
# Version: 1.6
# Function: 2D Monte Carlo

# Class - used to build isotopes information

import numpy as np

class Isotope:

    def __init__(self, name, atomic_mass, t_xs, e_xs, f_xs, nd):
        
        self.name = name

        self.mass = atomic_mass;

        self.alpha = ((atomic_mass - 1)/(atomic_mass+ 1)) ** 2;

        self.xs = np.zeros((3,3),float);
        self.ratio_react = np.zeros((3,3),float);
        self.macro_xs = np.zeros((3,1),float);
        self.scat_ratio = np.zeros((3,1),float);
        self.capt_ratio = np.zeros((3,1),float);
        self.fiss_ratio = np.zeros((3,1),float);
        self.absp_ratio = np.zeros((3,1),float);
        
        self.xs[0,] = t_xs;
        self.xs[1,] = e_xs;
        self.xs[2,] = f_xs;
        self.nd = nd; 
        
        for i in range(0,3):
            self.macro_xs[i] = np.sum(self.xs[i,]) * 1e-24 * self.nd;
            self.ratio_react[i,] = np.cumsum(self.xs[i,]) / np.sum(self.xs[i,]);

    def react_type(self, type):
        self.typee = type;
