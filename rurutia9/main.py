#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Author: Jiazhi Li
# Date: 12 Apr 2016
# Version: 1.4
# Function: 2D Monte Carlo
# Assumption:
#       No scattering - Only Capture(absorption and fission)
#       Isotopic outgoing angle
#       No reflection - Only leakage
#       Only transport mechanism no burnup
#       Two isotopes(H1 and U235)
#       Constant energy
#       ...

# Updates:
#       Implement sympy package
#       Implement periodic boundary conditions
#       Add scattering and absorption in H1
#       Add U238 and Pu239 isotopes

# test importing class file

import numpy as np
from test import Isotope

if __name__ == "__main__":
    a = Isotope('This is an isopotic name!');
    print a.name;
    b = Isotope('H');
    b.micro_xs([1,2,2,4],9.48*1e20);

    print b.capt_ratio;
