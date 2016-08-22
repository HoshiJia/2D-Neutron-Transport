#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Author: Jiazhi Li
# Date: 30 Apr 2016
# Version: 1.10
# Function: 2D Monte Carlo
# Assumption:
#       No scattering - Only Capture(absorption and fission)
#       Isotopic outgoing angle
#       No reflection - Only leakage
#       Only transport mechanism - no depletion
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
#       Record Fission Sites
#       Realize cycle iterations
#       Evaluate K effective value of each cycle and plot the trend
#       Abandon sympy package 
#       Expand the geometry to an arbitrary n x n pin-cells lattice

from __future__ import division
import sys
import numpy as np
from numpy import pi
import random
import math
##import matplotlib.pyplot as plt
import time

from isotope import Isotope
from neutron import Neutron
from geometry import Geometry
from copy import deepcopy
from matplotlib.backends.backend_pdf import PdfPages


def geo_intersc(neutron):


    # if position is inside the rectangle
    circle_flag = geo_circle(neutron.x_next,neutron.y_next);

    if neutron.energy >= 1 and neutron.energy <= 1e4 :
        g = 1; # resonance
    elif neutron.energy < 1:
        g = 0; # thermal
    else:
        g = 2; # faster
                        
    # Delta tracking
    if math.fabs(neutron.x_next) <= width and math.fabs(neutron.y_next) <= width:
        if circle_flag == 1:
            if random.random() < ratio_fuel_mod[g][0]:
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


def geo_circle(x1,y1,x0=0,y0=0,r=0.41):
    # check whether in the circle
    #c1 = Circle((x0,y0),r);
    #circle_in = c1.encloses_point((x1,y1));
    for i in range(tot_num):
        if (x1 - fuel.x[i]) ** 2 + (y1 - fuel.y[i]) ** 2 <= fuel.r[i] ** 2:
            circle_flag = 1;
            break;
        else:
            circle_flag = 2;
    return circle_flag;



def new_position(neutron):

    if neutron.energy >= 1 and neutron.energy <= 1e4 :
        g = 1; # resonance
    elif neutron.energy < 1:
        g = 0; # thermal
    else:
        g = 2; # faster

    # 0 - Delta tracking works; 1 - Delta tracking fails  
    travel_length = -math.log(random.random()) / mod_macro_xs[g][0];
    
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
        tan_t = math.sin(angle) / (1/neutron.iso_col.mass + math.cos(angle));
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

    return neutron;
    

def react_judge(neutron):

    while neutron.flag:

                #print neutron.x_rand, neutron.y_rand, neutron.x_next, neutron.y_next, neutron.angle;
                
                #p1 = Point(neutron.x_rand, neutron.y_rand);


                # (3)Pseudo Collision
                if neutron.flag == 2:
                    neutron.col = neutron.col + 1;

                    neutron.x_rand = neutron.x_next;
                    neutron.y_rand = neutron.y_next;
                            
                    neutron.iso_col='';
                

                # (2)exceed boundries
                elif neutron.flag == 3:
                    #neutron.col = neutron.col + 1;
                    while math.fabs(neutron.x_next) > width or math.fabs(neutron.y_next) > width:
                        if math.fabs(neutron.x_next) > math.fabs(neutron.y_next):
                            if neutron.x_next >= width:
                                neutron.x_next = neutron.x_next - 2 * width;
                                #neutron.y_next = neutron.y_next;
                            else:
                                neutron.x_next = neutron.x_next + 2 * width;
                                #neutron.y_next = neutron.y_next;
                        else:
                            if neutron.y_next >= b:
                                neutron.y_next = neutron.y_next - 2 * width;
                                #neutron.x_next = neutron.x_next;
                            else:
                                neutron.y_next = neutron.y_next + 2 * width;
                                #neutron.x_next = neutron.x_next;
                    neutron.iso_col='';
                    geo_intersc(neutron);
                    continue;
                        

                # (1) Accepted collision
                # check which reaction to happen (absorption or scattering)
                elif neutron.flag == 1:
                    neutron.col = neutron.col + 1;

                    if neutron.energy >= 1 and neutron.energy <= 1e4 :
                        g = 1; # resonance
                    elif neutron.energy < 1:
                        g = 0; # thermal
                    else:
                        g = 2; # faster
                                
                    # sample which isotope
                    rand_iso = random.random();
                    
                    # collision happens in the fuel;
                    if neutron.cur_domain == 1:
                        if rand_iso <= ratio_fuel_xs[g][0]:
                            # reacts with u235
                            neutron.iso_col = u235;
                        elif rand_iso > ratio_fuel_xs[g][1]:
                            # reacts with pu239
                            neutron.iso_col = pu239;
                        else:
                            # reacts with u238
                            neutron.iso_col = u238;

                    # scattering in moderator                            
                    elif neutron.cur_domain == 2:
                        if rand_iso <= ratio_mod_xs[g][0]:
                            # reacts with h1
                            neutron.iso_col = h1;
                        else:
                            # reacts with o16
                            neutron.iso_col = o16;
                            
                    else:
                        print ("Wrong Isotope Sampled!")
                            
                    rand_react = random.random();
                    
                    if rand_react <= neutron.iso_col.ratio_react[g][0]:
                        # fission
                        neutron.type = 2;
                        break; 

                    elif rand_react > neutron.iso_col.ratio_react[g][1]:
                        # scattering
                        neutron.type = 3;
                        neutron.x_rand = neutron.x_next;
                        neutron.y_rand = neutron.y_next;
                        neutron.angle = -1;
                    
                    else:
                        # absorption
                        neutron.type = 1;
                        break;


                if neutron.col>=1000:
                    print ('Lost tracking')
                    break; # avoid infinite loops

                
                neutron = new_position(neutron);
                geo_intersc(neutron);

    return neutron;



def main(neutron_cycle,cycle):

    
    # Fission and neutron source bank
    neutron_bank = np.zeros((cycle*neutron_cycle,8));
    # [No., site_x, site_y, travel_length, velocity_angle, reaction_type, number_neutron]

    rand_num = np.zeros((2,neutron_cycle),float);
    
    for n in range(neutron_cycle):
        rand_num[0][n] = random.uniform(-width,width);
        rand_num[1][n] = random.uniform(-width,width);

    rand_num_new = np.zeros((2,neutron_cycle),float);


    fission_num = neutron_cycle;
    
    for i in range(cycle):
        m = 0;
        #print rand_num;
        for j in range(neutron_cycle):

            row_num = neutron_cycle * i + j;
            neutron = Neutron(row_num,1);
            
            neutron_bank[row_num,1] = neutron.no;

            if j < fission_num:
                x_rand = rand_num[0][j];
                y_rand = rand_num[1][j];
            else:
                site_num = random.randint(0,fission_num-1);
                x_rand = rand_num[0][site_num];
                y_rand = rand_num[1][site_num];

            energy = 2*1e6; # unit: eV
            angle = -1;

            neutron_status = {'x_rand':x_rand,'y_rand':y_rand,'x_next':0,'y_next':0,'energy':1e6,'angle':-1,'iso_col':'','type':0};
            neutron.status(neutron_status);

            
            neutron_bank[row_num,2] = neutron.x_rand;
            neutron_bank[row_num,3] = neutron.y_rand;

            # print neutron_bank[row_num,2];
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
            # print neutron_bank[row_num,2];                
            neutron_bank[row_num,4] = neutron.x_next;
            neutron_bank[row_num,5] = neutron.y_next;
            neutron_bank[row_num,6] = neutron.col;
            neutron_bank[row_num,7] = neutron.type;

            if neutron.iso_col:
                x= neutron.iso_col.name;
            else:
                x= 'NaN';

            #print neutron.no, neutron.col, neutron_bank[row_num,2], neutron_bank[row_num,3], neutron_bank[row_num,4], neutron_bank[row_num,5], neutron.energy, x, neutron.type;
            # filter out fission sites
            if neutron.type == 2:
                rand_num_new[0][m] = deepcopy(neutron.x_next);
                rand_num_new[1][m] = deepcopy(neutron.y_next);
                #print rand_num_new[0][:]
                
                #try:
#                    print rand_num[0][:]
#                except:
#                    print 'xxx'
                m = m + 1;
                

            #print m;
            del neutron;
            #print rand_num;

        #print m;
        fission_num = m;
        fission_yield = 2.5;
        keff[i] = fission_num * fission_yield / neutron_cycle;
        print ('Keff = %f' % (fission_num * fission_yield / neutron_cycle));
        #del rand_num
        rand_num = deepcopy(rand_num_new);
        #del rand_num_new;

        
    end = time.time();
    print(end - start);


##    #pp = PdfPages('test.pdf')
##    
##    #colors = np.random.rand(neutron_cycle);
##    plt.figure(figsize=(8.4,8)); 
##    plt.subplot(2,2,1);
##    fig = plt.gcf();
##    for i in range(tot_num):
##        circle1 = plt.Circle((fuel.x[i],fuel.y[i]),fuel.r[i],alpha=0.1);
##        fig.gca().add_artist(circle1);
##    point1 = plt.scatter(neutron_bank[0:neutron_cycle,2], neutron_bank[0:neutron_cycle,3],alpha=0.5,lw = 0);
##    fig.gca().add_artist(point1);
##    plt.title('Neutron source - 1st cycle [start]');
##    plt.xlim([-width,width]);
##    plt.ylim([-width,width]);
##
##    plt.subplot(2,2,2);
##    fig = plt.gcf();
##    for i in range(tot_num):
##        circle1 = plt.Circle((fuel.x[i],fuel.y[i]),fuel.r[i],alpha=0.1);
##        fig.gca().add_artist(circle1);
##    c=range(4)
##    point1 = plt.scatter(neutron_bank[0:neutron_cycle,4], neutron_bank[0:neutron_cycle,5],c=neutron_bank[0:neutron_cycle,7],alpha=0.5,vmin=0, vmax=4,lw = 0);
##    fig.gca().add_artist(point1);
##    plt.title('Fission site - 1st cycle [end]');
##    plt.xlim([-width,width]);
##    plt.ylim([-width,width]);
##
##    plt.subplot(2,2,3);
##    fig = plt.gcf();
##    for i in range(tot_num):
##        circle1 = plt.Circle((fuel.x[i],fuel.y[i]),fuel.r[i],alpha=0.1);
##        fig.gca().add_artist(circle1);
##    point1 = plt.scatter(neutron_bank[neutron_cycle:neutron_cycle*2,2], neutron_bank[neutron_cycle:neutron_cycle*2,3],alpha=0.5,lw = 0);
##    fig.gca().add_artist(point1);
##    plt.title('Neutron source - 2nd cycle [start]');
##    plt.xlim([-width,width]);
##    plt.ylim([-width,width]);
##
##    plt.subplot(2,2,4);
##    fig = plt.gcf();
##    for i in range(tot_num):
##        circle1 = plt.Circle((fuel.x[i],fuel.y[i]),fuel.r[i],alpha=0.1);
##        fig.gca().add_artist(circle1);
##    point1 = plt.scatter(neutron_bank[neutron_cycle:neutron_cycle*2,4], neutron_bank[neutron_cycle:neutron_cycle*2,5],c=neutron_bank[neutron_cycle:neutron_cycle*2,7],alpha=0.5,vmin=0, vmax=4,lw = 0,s=20);
##    fig.gca().add_artist(point1);
##    plt.title('Fission site - 2nd cycle [end]');
##    plt.xlim([-width,width]);
##    plt.ylim([-width,width]);
##
##    #pp.savefig()
##    #pp.close()
##    
##    #plt.figure(2);
##    #plt.plot(keff);
##    
##    
##    plt.show();

    #time.sleep(1);



if __name__ == "__main__":
    start = time.time()
    neutron_cycle = 5000;
    cycle = 10;

    r = 0.41;
    b = 0.627;

    pin_num = 3;
    tot_num = pin_num ** 2;

    width = b * pin_num;
    fuel = Geometry('circle');

    fuel.circle(r, 2*b, pin_num);
    keff = np.zeros((cycle,1),float);
    en_gp = 3;
    
##    h1 = Isotope('h1', 1, [0,0.2,20], [0,0.00004,4],6.7358*1e22); # calculate macro cross-section
##    o16 = Isotope('o16', 16, [0,0.0001,4], [0,0.00000003,3],3.3679*1e22);
##    u235 = Isotope('u235', 235, [530,99,10], [1,0.09,4],9.3472*1e20);
##    u238 = Isotope('u238', 238, [0.00002,2,9], [0.3,0.07,5],2.1523*1e22);
##    pu239 = Isotope('pu239', 239, [748,269,8], [2,0.05,5],1*1e18);
    h1 = Isotope('h1', 1, [0,0.294,28.966],[0,0.149,240.23], [0,0.00004,4],6.7358*1e22); # calculate macro cross-section
    o16 = Isotope('o16', 16, [0,0.0001,4], [0,0.0001,4],[0,0.00000003,3],3.3679*1e22);
    u235 = Isotope('u235', 235, [504.8,86.7,15.98], [271.5,131.97,152.82],[1.219,0.095,6.3],9.3472*1e20);
    u238 = Isotope('u238', 238, [0.00001,2.4,9.4], [0.002,277.7,319.06],[0.3,0.07,4.825],2.1523*1e22);
    pu239 = Isotope('pu239', 239, [699.34,274.32,7.90], [289.36,184.06,155.87], [1.369,1.800,0.065],0);

    
    macro_xs = np.zeros((en_gp,5),float);
    mod_macro_xs = np.zeros((en_gp,1),float);
    fuel_macro_xs = np.zeros((en_gp,1),float);
    ratio_mod_xs = np.zeros((en_gp,2),float);
    ratio_fuel_xs = np.zeros((en_gp,3),float);
    ratio_fuel_mod = np.zeros((en_gp,1),float);
    for i in range(0,en_gp):
        macro_xs[i,] = [h1.macro_xs[i], o16.macro_xs[i], u235.macro_xs[i], u238.macro_xs[i], pu239.macro_xs[i]];
        mod_macro_xs[i] = np.sum(macro_xs[i][0:2]);
        fuel_macro_xs[i] = np.sum(macro_xs[i][2:5]);
        ratio_mod_xs[i,] = np.cumsum(macro_xs[i][0:2]) / mod_macro_xs[i];
        ratio_fuel_xs[i,] = np.cumsum(macro_xs[i][2:5]) / fuel_macro_xs[i];
        ratio_fuel_mod[i] = fuel_macro_xs[i] / mod_macro_xs[i];

    #print mod_macro_xs,fuel_macro_xs;
    
    main(neutron_cycle,cycle);
    
