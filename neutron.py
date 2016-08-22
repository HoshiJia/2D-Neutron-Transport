#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Author: Jiazhi Li
# Date: 15 Apr 2016
# Version: 1.6
# Function: 2D Monte Carlo

# Class - used to build neutron information


class Neutron:

    def __init__(self, no, n):
        
        self.no = no;
        self.col = n;

    def status(self, neutron_status):

        self.x_rand = neutron_status['x_rand'];
        self.y_rand = neutron_status['y_rand'];

        self.x_next = neutron_status['x_next'];
        self.y_next = neutron_status['y_next'];

        self.angle = neutron_status['angle'];
        self.energy = neutron_status['energy'];
        self.iso_col = neutron_status['iso_col'];
        self.type = neutron_status['type'];

        self.x_static = neutron_status['x_rand'];
        self.y_static = neutron_status['y_rand'];

    def intersect(self, inter_sc, flag):
        self.flag = flag;
        self.pos = inter_sc;

    def domain(self, domain):
        self.cur_domain = domain;

    def sample(self, resample):
        self.resample = resample;

    
