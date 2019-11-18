########################################
# Usage: 'python Prep.py ver prefix'
#
# All simulation scripts can be found by the names Prep.py, Min.py, Equ.py and Prod.py. You can see the protocol by reading them.

########################################
# Importing python libraries

import re
import sys
import os
import fileinput

########################################
# Defining variables
ver = sys.argv[1]       # mpi parameters
prefix = sys.argv[2]    # molecule name (e.g. ben)

###########################################
# create input coordinates and topology
if not os.path.exists(prefix + ".gro"):
	pdb2gmxcommand = "pdb2gmx_" + ver + " -ff argromos -water spc -f " + prefix + ".pdb -p " + prefix + ".top -o " + prefix + ".gro"
	os.system(pdb2gmxcommand)
else:
	print "\n>>> Topology created. Moving on..."

###########################################
# setup periodic box
if not os.path.exists(prefix + ".pbc.gro"):
	setupcommand = "editconf_" + ver + " -f " + prefix + ".gro -c -o " + prefix + ".pbc.gro -bt cubic -box 2"
	os.system(setupcommand)
else:
	print "\n>>> Box created. Moving on..."

###########################################
# Copying files for Gas-Phase simulation
if not os.path.exists(prefix + ".1.top"):
	os.system("cp " + prefix + ".top " + prefix + ".1.top")

	setupcommand = "editconf_" + ver + " -f " + prefix + ".gro -c -o " + prefix + ".1.pbc.gro -bt cubic -box 4"
	os.system(setupcommand)
else:
	print "\n Gas-phase input created. Moving on..."
