########################################
# Usage: 'python LIQ-sim.py ver prefix'
#
# This script was written to run Liquid phase simulations. This protocol was extracted from "Force Field Benchmark of Organic Liquids: Density, Enthalpy of Vaporization, Heat Capacities, Surface Tension, Isothermal Compressibility, Volumetric Expansion Coefficient, and Dielectric Constant (Caleman, Carl. 2011);
#
# All simulation scripts can be found by the names Prep.py, Min.py, Equ.py and Prod.py. You can see the protocol by reading them.

########################################
# Importing python libraries
import os
import sys
import datetime
########################################
# Defining variables

ver = sys.argv[1]
liq = sys.argv[2]

########################################
# Submitting Liquid system simulation

time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print ">>> Starting time = " + time

os.system("python Prep.py " + ver + " " + liq)
os.system("python Min.py " + ver + " " + liq)
os.system("python Equ.py " + ver + " " + liq)
os.system("python 3press.py " + ver + " " + liq)
os.system("python 3temp.py " + ver + " " + liq)
os.system("python LIQ.py " + ver + " " + liq)
os.system("python SD.py " + ver + " " + liq)

print ">>> Ending time = " + time
print "Your simulation has finished. Congratulations!"
