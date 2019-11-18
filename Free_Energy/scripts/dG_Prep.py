########################################
# Usage: 'python Prep.py ver prefix'
#
# All simulation scripts can be found by the names Prep.py, Min.py, Equ.py and Prod.py. You can see the protocol by reading them.

########################################
# Importing python libraries

import sys
import os

########################################
# Defining variables
ver = sys.argv[1]       # mpi parameters
prefix = sys.argv[2]    # molecule name (e.g. ben)

###########################################
# create input coordinates and topology
if not os.path.exists(prefix + ".gro"):
	pdb2gmxcommand = "gmx_" + ver + " pdb2gmx -ff argromos -water spc -f " + prefix + ".pdb -p " + prefix + ".top -o " + prefix + ".gro"
	os.system(pdb2gmxcommand)
else:
	print "\n>>> Topology created. Moving on..."

###########################################
# setup periodic box
if not os.path.exists(prefix + ".pbc.gro"):
	setupcommand = "gmx_" + ver + " editconf -f " + prefix + ".gro -c -o " + prefix + ".pbc.gro -bt cubic -box 3"
	os.system(setupcommand)
else:
	print "\n>>> Box created. Moving on..."

###########################################
# Inserting solvent
if not os.path.exists(prefix + ".sol0.gro"):
	setupcommand = "gmx_" + ver + " solvate -cp " + prefix + ".pbc.gro -cs spc216.gro -o " + prefix + ".sol0.gro -p " + prefix + ".top"
	os.system(setupcommand)
else:
	print "\n>>> Solvent inserted. Moving on..."

###########################################
# setup periodic box
if not os.path.exists(prefix + ".sol.gro"):
	setupcommand = "gmx_" + ver + " editconf -f " + prefix + ".sol0.gro -c -o " + prefix + ".sol.gro -bt cubic -density 997"
	os.system(setupcommand)
else:
	print "\n>>> Density fixed. Moving on..."


def Minimization0():
    inmin0 = "in.min0.mdp"
    inmin0_h = open(inmin0, 'w')
    content = """; Run control
    integrator               = steep
    nsteps                   = 1000
    ; EM criteria and other stuff
    emtol                    = 100
    emstep                   = 0.01
    niter                    = 20
    nbfgscorr                = 10
    ; Output control
    nstlog                   = 1
    nstenergy                = 1
    ; Neighborsearching and short-range nonbonded interactions
    cutoff-scheme            = verlet
    nstlist                  = 1
    ns_type                  = grid
    pbc                      = xyz
    rlist                    = 1.2
    ; Electrostatics
    coulombtype              = PME
    rcoulomb                 = 1.2
    ; van der Waals
    vdwtype                  = cutoff
    vdw-modifier             = potential-switch
    rvdw-switch              = 1.0
    rvdw                     = 1.2
    ; Apply long range dispersion corrections for Energy and Pressure
    DispCorr                  = EnerPres
    ; Spacing for the PME/PPPM FFT grid
    fourierspacing           = 0.12
    ; EWALD/PME/PPPM parameters
    pme_order                = 6
    ewald_rtol               = 1e-06
    epsilon_surface          = 0
    ; Temperature and pressure coupling are off during EM
    tcoupl                   = no
    pcoupl                   = no
    ; No velocities during EM
    gen_vel                  = no
    ; options for bonds
    constraints              = h-bonds  ; we only have C-H bonds here
    ; Type of constraint algorithm
    constraint-algorithm     = lincs
    ; Do not constrain the starting configuration
    continuation             = no
    ; Highest order in the expansion of the constraint coupling matrix
    lincs-order              = 12

        """
    inmin0_h.write(content)
    inmin0_h.close()

    if not os.path.exists(prefix + ".min0.gro"):
    	min1gromppcommand = "gmx_" + ver + " grompp -f in.min0.mdp -c " + prefix + ".sol.gro -p " + prefix + ".top -o " + prefix + ".min0.tpr -po " + prefix + ".min1.mdout.mdp -maxwarn 1"
    	os.system(min1gromppcommand)

    	min1command = "gmx_" + ver + " mdrun -deffnm " + prefix + ".min0"
    	os.system(min1command)
    else:
    	print "\n>>> Minimization 0 completed!\n"

Minimization0()
