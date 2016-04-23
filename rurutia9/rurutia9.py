#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Author: Jiazhi Li
# Date: 21 Apr 2016
# Version: 1.8
# Function: 2D Monte Carlo
# Assumption:
#       No scattering - Only Capture(absorption and fission)
#       Isotopic outgoing angle
#       No reflection - Only leakage
#       Only transport mechanism - no burnup
#       Two isotopes(H1 and U235)
#       Constant energy
#       ...

# Updates:
#       Implement sympy package
#       Implement periodic boundary conditions
#       Add scattering and absorption
#       Add U238 and Pu239 isotopes
#       Two energy groups and energy loss due to scattering
#       Optimize periodic algorithm
#       Implement Delta Tracking method

from __future__ import division
import sys
import numpy as np
from numpy import pi
import random
import math
from sympy import *
import matplotlib.pyplot as plt
import time

from isotope import Isotope
from neutron import Neutron



def geo_intersc(neutron):
    
    
    p1, p2 = Point(neutron.x_rand,neutron.y_rand), Point(neutron.x_next,neutron.y_next);
    
    r1 = Polygon((-b,b),(-b,-b),(b,-b),(b,b));

    # if position is inside the rectangle
    circle_flag = geo_circle(p2.x,p2.y);

    # Delta tracking
    if r1.encloses_point(p2):
        if circle_flag == 1:
            if random.random() < ratio_fuel_mod:
                neutron.flag = 1;
                neutron.domain(1);
                # accepted
            else:
                neutron.flag = 2;
                neutron.domain(1);
        else:
            neutron.flag = 1;
            neutron.domain(2);

    else:
        neutron.flag = 3;
        neutron.domain(3);

    # flag -1. Accepted collision -2. Pseudo collision  -3. Exceed the boundary

        
    return neutron;


def geo_circle(x1,y1,x0=0,y0=0,r=0.45):
    # check whether in the circle
    c1 = Circle((x0,y0),r);
    circle_in = c1.encloses_point((x1,y1));
    if circle_in:
        circle_flag = 1;
    else:
        circle_flag = 2;
    return circle_flag;



def new_position(neutron):

    # 0 - Delta tracking works; 1 - Delta tracking fails  
    travel_length = -math.log(random.random()) / mod_macro_xs;
    travel_length = float("{0:.6f}".format(travel_length));


    # generate random angle
    if neutron.angle == -1 : # incident angle or inelastic collision angle in COM
        #angle = N(2 * pi * random.random(),6);
        angle = 2 * pi * random.random();
    else: # crossing boundary
        #angle = N(neutron.angle,6);
        angle = neutron.angle;
        
    if neutron.iso_col == '': # -1: initial neutron
        u = math.cos(angle); # incident or collision
        v = math.sin(angle);
    else:
        u_c = math.cos(angle); # collision angle in COM [0, 2*pi]
        tan_t = math.sin(angle) / (1/neutron.iso_col.no + math.cos(angle));
        angle_real = math.atan(tan_t);
        #u_l = (1 + neutron.iso_col.no * u_c) / math.sqrt(1 + 2 * neutron.iso_col.no * u_c + neutron.iso_col.no ** 2);
        neutron.energy = neutron.energy * 1/2 * ((1 + neutron.iso_col.alpha) + (1 - neutron.iso_col.alpha) * u_c);

        u = math.cos(angle);
        v = math.sin(angle);
            
    neutron.iso_col == '';
    neutron.angle = angle;
    
    # fission position
    
    neutron.x_next = neutron.x_rand + travel_length * u;
    neutron.y_next = neutron.y_rand + travel_length * v;

    neutron.x_next = float("{0:.6f}".format(neutron.x_next));
    neutron.y_next = float("{0:.6f}".format(neutron.y_next));
    #neutron.x_next = "%.3f" % neutron.x_next;
    #neutron.y_next = "%.3f" % neutron.y_next;    

    return neutron;
    

def react_judge(neutron):

    while neutron.flag:

                #print neutron.x_rand, neutron.y_rand, neutron.x_next, neutron.y_next, neutron.angle;
                
                p1 = Point(neutron.x_rand, neutron.y_rand);


                # (3)Pseudo Collision
                if neutron.flag == 2:
                    neutron.col = neutron.col + 1;

                    neutron.x_rand = neutron.x_next;
                    neutron.y_rand = neutron.y_next;
                            
                    neutron.iso_col='';
                

                # (2)exceed boundries
                elif neutron.flag == 3:
                    #neutron.col = neutron.col + 1;
                    while math.fabs(neutron.x_next) > b or math.fabs(neutron.y_next) > b:
                        if math.fabs(neutron.x_next) > math.fabs(neutron.y_next):
                            if neutron.x_next >= b:
                                neutron.x_next = neutron.x_next - 2 * b;
                                #neutron.y_next = neutron.y_next;
                            else:
                                neutron.x_next = neutron.x_next + 2 * b;
                                #neutron.y_next = neutron.y_next;
                        else:
                            if neutron.y_next >= b:
                                neutron.y_next = neutron.y_next - 2 * b;
                                #neutron.x_next = neutron.x_next;
                            else:
                                neutron.y_next = neutron.y_next + 2 * b;
                                #neutron.x_next = neutron.x_next;
                    neutron.iso_col='';
                    geo_intersc(neutron);
                    continue;
                        

                # (1) Accepted collision
                # check which reaction to happen (absorption or scattering)
                elif neutron.flag == 1:
                    neutron.col = neutron.col + 1;
                    
                    # collision happens in the fuel;
                    if neutron.cur_domain == 1:

                        # scattering in fuel
                        # sample which isotope
                        samp = random.random();
                        if samp <= ratio_fuel_xs[0]:
                            # reacts with u235
                            if random.random() > u235.scat_ratio:
                                neutron.iso_col = u235;
                                
                                break; # only way to jump loop - absorption
                                # continued to supplement
                                
                            else:
                                # scattering
                                neutron.x_rand = neutron.x_next;
                                neutron.y_rand = neutron.y_next;
                                neutron.angle = -1;
                                neutron.iso_col = u235;
                                
                        elif samp > ratio_fuel_xs[1]:
                            # reacts with pu239
                            if random.random() > pu239.scat_ratio:
                                neutron.iso_col = pu239;
                                break; # only way to jump loop - absorption
                                # continued to supplement 
                            else:
                                # scattering
                                neutron.x_rand = neutron.x_next;
                                neutron.y_rand = neutron.y_next;
                                neutron.angle = -1;
                                neutron.iso_col = pu239;
                        else:
                            # reacts with u238
                            if random.random() > u238.scat_ratio:
                                neutron.iso_col = u238;
                                break; # only way to jump loop - absorption
                                # continued to supplement 
                            else:
                                # scattering
                                neutron.x_rand = neutron.x_next;
                                neutron.y_rand = neutron.y_next;
                                neutron.angle = -1;
                                neutron.iso_col = u238;
                        #print neutron.iso_col;
                    elif neutron.cur_domain == 2:

                        # scattering in moderator
                        # sample from which isotope
                        samp = random.random();
                        if samp <= ratio_mod_xs[0]:
                            # reacts with H1
                            if random.random() > h1.scat_ratio:
                                break; # only way to jump loop - absorption
                            else:
                                #neutron.x_rand = N(neutron.x_next,6);
                                #neutron.y_rand = N(neutron.y_next,6);
                                neutron.x_rand = neutron.x_next;
                                neutron.y_rand = neutron.y_next;
                                neutron.angle = -1;
                                neutron.iso_col = h1;
                        else:
                            # reacts with o16
                            if random.random() > o16.scat_ratio:
                                break; # only way to jump loop - absorption
                            else:
                                #neutron.x_rand = N(neutron.x_next,6);
                                #neutron.y_rand = N(neutron.y_next,6);
                                neutron.x_rand = neutron.x_next;
                                neutron.y_rand = neutron.y_next;
                                neutron.angle = -1;
                                neutron.iso_col = o16;
                    else:
                        print "Wrong Isotope Sampled!"
                    
                if neutron.col>=10:
                    print 'Lost tracking'
                    break; # avoid infinite loops

                
                neutron = new_position(neutron);
                geo_intersc(neutron);

    return neutron;



def main(neutron_cycle,cycle):

    
    # Fission and neutron source bank
    neutron_bank = np.zeros((neutron_cycle,7));
    # [No., site_x, site_y, travel_length, velocity_angle, reaction_type, number_neutron]
    
    
    for i in range(cycle):
        for j in range(neutron_cycle):
            
            neutron = Neutron(j,1);
            
            neutron_bank[j,1] = neutron.no;

            x_rand = 2*b*random.random()-b;
            y_rand = 2*b*random.random()-b;
            #x_rand = "%.6f" % rand_num_1;
            #y_rand = "%.6f" % rand_num_2;

            #x_rand = float(x_rand);
            #y_rand = float(y_rand);

            energy = 1e6; # unit: eV
            angle = -1;

            neutron_status = {'x_rand':x_rand,'y_rand':y_rand,'x_next':0,'y_next':0,'energy':1e6,'angle':-1,'iso_col':''};
            neutron.status(neutron_status);


            neutron_bank[j,2] = neutron.x_rand;
            neutron_bank[j,3] = neutron.y_rand;


            # generate the first collision position
            circle_flag = geo_circle(neutron.x_rand, neutron.y_rand);

            if circle_flag == 1:
                neutron.domain(1);
            else:
                neutron.domain(2);
            
            neutron = new_position(neutron);
            
            geo_intersc(neutron);
            
            
            #neutron_bank[j,4] = neutron['x_next'];
            #neutron_bank[j,5] = neutron['y_next'];
            
            # loop core
            react_judge(neutron);
                
            neutron_bank[j,4] = neutron.x_next;
            neutron_bank[j,5] = neutron.y_next;
            neutron_bank[j,6] = neutron.col;

            #if neutron.iso_col:
            #    x= neutron.iso_col.name;
            #else:
            #    x= 'NaN'

            #print neutron.no, neutron.col, neutron.x_rand, neutron.y_rand, neutron.x_next, neutron.y_next, neutron.energy, x;
            print neutron.no;
            del neutron;
            
            #print n;
            #print j, n;
    end = time.time();
    print(end - start);
    plt.subplot(1,2,1);
    circle1 = plt.Circle((0,0),r,alpha=0.1);
    point1 = plt.scatter(neutron_bank[:,2], neutron_bank[:,3]);
    fig = plt.gcf();
    fig.gca().add_artist(circle1);
    fig.gca().add_artist(point1);
    plt.title('Neutron source');
    plt.xlim([-b,b]);
    plt.ylim([-b,b]);

    plt.subplot(1,2,2);
    circle1 = plt.Circle((0,0),r,alpha=0.1);
    point1 = plt.scatter(neutron_bank[:,4], neutron_bank[:,5]);
    fig = plt.gcf();
    fig.gca().add_artist(circle1);
    fig.gca().add_artist(point1);
    plt.title('Fission site');
    plt.xlim([-b,b]);
    plt.ylim([-b,b]);
    
    plt.show();
    #time.sleep(1);



if __name__ == "__main__":
    start = time.time()
    neutron_cycle = 2000;
    cycle = 1;

    r = 0.45;
    b = 0.627;
    
    h1 = Isotope('h1',[1,0,0.2,20],6.7358*1e22); # calculate macro cross-section
    o16 = Isotope('o16',[16,0,0.0001,4],3.3679*1e22);
    u235 = Isotope('u235',[235,530,99,10],9.3472*1e20);
    u238 = Isotope('u238',[238,0.00002,2,9],2.1523*1e22);
    pu239 = Isotope('pu239',[239,748,269,8],1*1e18);
    
    macro_xs = [h1.macro_xs, o16.macro_xs, u235.macro_xs, u238.macro_xs, pu239.macro_xs];
    mod_macro_xs = np.sum(macro_xs[0:2]);
    fuel_macro_xs = np.sum(macro_xs[2:5]);
    ratio_mod_xs = np.cumsum(macro_xs[0:2]) / mod_macro_xs;
    ratio_fuel_xs = np.cumsum(macro_xs[2:5]) / fuel_macro_xs;
    ratio_fuel_mod = fuel_macro_xs / fuel_macro_xs;
    
    main(neutron_cycle,cycle);
    
