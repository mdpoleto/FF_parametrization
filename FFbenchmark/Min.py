########################################
# Usage: 'python Min.py ver prefix'
#
# All simulation scripts can be found by the names Prep.py, Min.py, Equ.py, Prod.py and SD.py. You can see the protocol by reading them.

########################################
# Importing python libraries
import re
import sys
import os

########################################
# Defining variables
ver = sys.argv[1]       # mpi parameters
prefix = sys.argv[2]    # molecule name (e.g. ben)

#########################################
def change_line(file,param):

	f_in = open(file,"r")
	f_out = open(file[:-4]+".out.top","w")
	line = " dad w "
	new_line = ""
	while True:
		line = f_in.readline()
		if not line: break
		
		if len(line.split()) > 0:
			if line.split()[0] == "SOL":
				
				if line.split()[1] == "1":
					new_line = "SOL                 "+str(param)
					
					f_out.write(new_line)
				else:
					f_out.write(line)
			else:
				f_out.write(line)
		else:
			f_out.write(line)
	f_in.close()
	f_out.close()

##########################################
# Stacking boxes for Liquid-Phase simulation 

if not os.path.exists(prefix + ".stack.gro"):
	genconfcommand = "genconf_" + ver + " -f " + prefix + ".pbc.gro -o " + prefix + ".stack.gro -nbox 5 5 5"
	os.system(genconfcommand)

	change_line(prefix + ".top", 125)
	os.system("rm " + prefix + ".top")
	os.system("mv " + prefix + ".out.top " + prefix + ".top")
else:
	print "\n>>> Box stacked succesfully. Moving on..."

############################################
# solvated system energy minimization 1

inmin1 = "in.min1.mdp"
inmin1_h = open(inmin1, 'w')
content = """; Define can be used to control processes
    
    define          = -DFLEXIBLE
    
    ; Parameters describing what to do, when to stop and what to save
    integrator      = steep
    emtol           = 100.0
    emstep          = 0.01
    nsteps          = 10000
    nstenergy       = 1
    energygrps      = System
    
    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
    nstlist         = 1
    ns_type         = grid
    coulombtype     = PME
    rlist           = 1.0
    rcoulomb        = 1.0
    rvdw            = 1.0
    constraints     = none
    pbc             = xyz
    """
inmin1_h.write(content)
inmin1_h.close()

if not os.path.exists(prefix + ".min1.gro"):
	min1gromppcommand = "grompp_" + ver + " -f in.min1.mdp -c " + prefix + ".stack.gro -p " + prefix + ".top -o " + prefix + ".min1.tpr -po " + prefix + ".min1.mdout.mdp -maxwarn 1"
	os.system(min1gromppcommand)

	min1command = "mdrun_" + ver + " -deffnm " + prefix + ".min1"
	os.system(min1command)
else:
	print "\n>>> Minimization 1 completed!\n"
############################################
# solvated system energy minimization 2

inmin2 = "in.min2.mdp"
inmin2_h = open(inmin2, 'w')
content = """; Define can be used to control processes
    
    define          = -DFLEXIBLE
    
    ; Parameters describing what to do, when to stop and what to save
    integrator      = cg
    emtol           = 100.0
    emstep          = 0.01
    nsteps          = 10000
    nstenergy       = 1
    energygrps      = System
    
    ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
    nstlist         = 1
    ns_type         = grid
    coulombtype     = PME
    rlist           = 1.0
    rcoulomb        = 1.0
    rvdw            = 1.0
    constraints     = none
    pbc             = xyz
    """
inmin2_h.write(content)
inmin2_h.close()

if not os.path.exists(prefix + ".min2.gro"):
	min2gromppcommand = "grompp_" + ver + " -f in.min2.mdp -c " + prefix + ".min1.gro -p " + prefix + ".top -o " + prefix + ".min2.tpr -po " + prefix + ".min2.mdout.mdp"
	os.system(min2gromppcommand)

	min2command = "mdrun_" + ver + " -deffnm " + prefix + ".min2"
	os.system(min2command)
else:
	print "\n>>> Minimization 2 completed!\n"
