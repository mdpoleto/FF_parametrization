########################################
# Usage: 'python SD.py ver prefix'
#
# All simulation scripts can be found by the names Prep.py, Min.py, Equ.py, Prod.py and SD.py. You can see the protocol by reading them.

########################################
# Importing python libraries
import re
import sys
import os

########################################
# Defining variables

ver = sys.argv[1]       # gromacs version (e.g. 504)
prefix = sys.argv[2]    # molecule name (e.g. ben)

########################################
# Grepping properties

db_directory = '/home/marcelodepolo/Heterocycles/smGROMOS-prop.dab'

def grep_prop(file,prefix):
	f_in = open(file,"r")
	line = ""
	
	while True:
		line = f_in.readline()
		if not line: break
		
		if len(line.split()) > 0:
			if line.split()[1] == prefix.upper(): 
				name = line.split()[0]
				code = line.split()[1]
				T = line.split()[2]
				DF = line.split()[3]
				ExD = line.split()[4]
				ExV = line.split()[5]
				ExE = line.split()[6]
				ExC = line.split()[7]
			else:
				pass
		else:
			pass
	return name, code, T, DF, ExD, ExV, ExC, ExE

prop = grep_prop(db_directory,prefix)

name = str(prop[0])
code = str(prop[1])
T = float(prop[2])   # K
DF = int(prop[3])  
if prop[6] == 'nan':
	ExC = 0.5            # When there is no data for Isoth. Compressibility, default = 0.5 1/GPa
else:
	ExC = float(prop[6]) # 1GPa = 9.5e-5 bar^-1
	
comp = float(ExC/(10**4))

def mdp_equ(ver,prefix):
	############################################
	# simulation MD, production run

	inmd = "in.SDequ.mdp"
	inmd_h = open(inmd, 'w')
	content = """; RUN CONTROL PARAMETERS
integrator               = sd
dt                       = 0.002
nsteps                   = 2500000
comm_mode                = none

; OUTPUT CONTROL OPTIONS
nstxout         = 0
nstvout         = 0
nstenergy       = 1000
nstlog          = 1000
nstxtcout       = 0
xtc-precision   = 1000
energygrps      = System

; NEIGHBORSEARCHING PARAMETERS
nstlist         = 10
ns-type         = Grid
pbc             = xyz
rlist           = 1.1
cutoff-scheme   = Verlet

; OPTIONS FOR ELECTROSTATICS AND VDW
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.12
rcoulomb        = 1.1
vdw-type        = PME
rvdw            = 1.1

; Temperature coupling
Tcoupl          = nose-hoover
tc-grps         = System
tau_t           = 0.1
ref_t           = """

	temp = str(T)

	content2 = """
; Dispersion correction
DispCorr                 = EnerPres ; account for vdw cut-off

; Pressure coupling
Pcoupl                   = no	

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                 = yes    ; Assign velocities to particles by taking them randomly from a Maxwell distribution
gen_temp                = """

	content3 = """    ; Temperature to generate corresponding Maxwell distribution
gen_seed                = -1     ; Seed for (semi) random number generation. Different numbers give different sets of velocities

; OPTIONS FOR BONDS
constraints             = all-bonds
continuation            = no
constraint_algorithm    = lincs
lincs_iter              = 2
lincs_order             = 4
	"""
	inmd_h.write(content)
	inmd_h.write(temp)
	inmd_h.write(content2)
	inmd_h.write(temp)
	inmd_h.write(content3)
	inmd_h.close()

def mdp_prod(ver,prefix):
	############################################
	# simulation MD, production run

	inmd = "in.SD.mdp"
	inmd_h = open(inmd, 'w')
	content = """; RUN CONTROL PARAMETERS
integrator               = sd
dt                       = 0.002
nsteps                   = 50000000
comm_mode                = none

; OUTPUT CONTROL OPTIONS
nstxout         = 0
nstvout         = 0
nstenergy       = 1000
nstlog          = 1000
nstxtcout       = 0
xtc-precision   = 1000
energygrps      = System

; NEIGHBORSEARCHING PARAMETERS
nstlist         = 10
ns-type         = Grid
pbc             = xyz
rlist           = 1.1
cutoff-scheme   = Verlet

; OPTIONS FOR ELECTROSTATICS AND VDW
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.12
rcoulomb        = 1.1
vdw-type        = PME
rvdw            = 1.1

; Temperature coupling
Tcoupl          = nose-hoover
tc-grps         = System
tau_t           = 0.1
ref_t           = """

	temp = str(T)

	content2 = """
; Dispersion correction
DispCorr                 = EnerPres ; account for vdw cut-off

; Pressure coupling
Pcoupl                   = no	

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                 = yes    ; Assign velocities to particles by taking them randomly from a Maxwell distribution
gen_temp                = """

	content3 = """    ; Temperature to generate corresponding Maxwell distribution
gen_seed                = -1     ; Seed for (semi) random number generation. Different numbers give different sets of velocities

; OPTIONS FOR BONDS
constraints             = all-bonds
continuation            = no
constraint_algorithm    = lincs
lincs_iter              = 2
lincs_order             = 4
	"""
	inmd_h.write(content)
	inmd_h.write(temp)
	inmd_h.write(content2)
	inmd_h.write(temp)
	inmd_h.write(content3)
	inmd_h.close()

if not os.path.exists(prefix + ".SDequ.cpt"):
	
	mdp_equ(ver,prefix)
	
	gromppcommand = "grompp_" + ver + " -f in.SDequ.mdp -c " + prefix + ".1.pbc.gro -p " + prefix + ".1.top -o " + prefix + ".SDequ.tpr -po " + prefix + ".SDequ.mdout.mdp -maxwarn 1"
	os.system(gromppcommand)

	mdruncommand = "mdrun_" + ver + " -deffnm " + prefix + ".SDequ"
	os.system(mdruncommand)
else:
   mdruncommand = "mdrun_" + ver + " -deffnm " + prefix + ".SDequ -cpi " + prefix + ".SDequ.cpt -append"
   os.system(mdruncommand)

if not os.path.exists(prefix + ".SD.cpt"):
	
	mdp_prod(ver,prefix)
	
	gromppcommand = "grompp_" + ver + " -f in.SD.mdp -c " + prefix + ".SDequ.gro -p " + prefix + ".1.top -o " + prefix + ".SD.tpr -po " + prefix + ".SD.mdout.mdp -maxwarn 1"
	os.system(gromppcommand)

	mdruncommand = "mdrun_" + ver + " -deffnm " + prefix + ".SD"
	os.system(mdruncommand)
else:
   mdruncommand = "mdrun_" + ver + " -deffnm " + prefix + ".SD -cpi " + prefix + ".SD.cpt -append"
   os.system(mdruncommand)

