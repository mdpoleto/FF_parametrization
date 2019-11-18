########################################
# Usage: 'python Prod.py ver prefix'
#
# All simulation scripts can be found by the names Prep.py, Min.py, Equ.py, Prod.py and SD.py. You can see the protocol by reading them.

########################################
# Importing python libraries
import re
import sys
import os

########################################
# Setting Sys Arguments
ver = sys.argv[1]       # version
prefix = sys.argv[2]    # molecule name (e.g. ben)

#######################################
# Constant values
Kb = 0.0083 # kJ/mol/K
N = 1000 # Number of molecules in final system
threshold = 0.5 # J/mol/ns/DegFree (Total-Energy drift)

###############################################
# Reading Molecules Libraries (DF = Degree of freedom; ExD = Experimental Density in g/cm^3; ExV = Experimental Vaporization Enthalpy kJ/mol)

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
				ExDc = line.split()[8]
				ExDm = line.split()[9]
				ExCp = line.split()[10]
				ExCv = line.split()[11]
				MW = line.split()[12]
			else:
				pass
		else:
			pass
	return name, code, T, DF, ExD, ExV, ExE, ExC, ExDc, ExDm, ExCp, ExCv, MW

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

ExDc = float(prop[8])
diec_cutlow = 5
diec_cuthigh = 60

if ExDc <= diec_cutlow:
        dielectric = 5
elif ExDc > diec_cutlow and ExDc < diec_cuthigh:
        dielectric = ExDc
if ExDc >= diec_cutlow:
        dielectric = diec_cuthigh

########################################
# Defining Function

def create_mdp(ver,prefix):
	############################################
	# simulation mdp

	inmd = "in.LIQ.mdp"
	inmd_h = open(inmd, 'w')
	content = """; RUN CONTROL PARAMETERS
integrator               = md
tinit                    = 0
dt                       = 0.002
nsteps                   = 50000000
; For exact run continuation or redoing part of a run
init_step                = 0
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 1

; LANGEVIN DYNAMICS OPTIONS
; Friction coefficient (amu/ps) and random seed
bd-fric                  = 0
ld-seed                  = 1993

; ENERGY MINIMIZATION OPTIONS
; Force tolerance and initial step-size
emtol                    = 0.1
emstep                   = 0.01
; Max number of iterations in relax_shells
niter                    = 20
; Step size (ps^2) for minimization of flexible constraints
fcstep                   = 0
; Frequency of steepest descents steps when doing CG
nstcgsteep               = 1000
nbfgscorr                = 10

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxtcout                = 1000
; Output frequency for energies to log file and energy file
nstlog                   = 100
nstenergy                = 100
; Output frequency and precision for xtc file
xtc-precision            = 1000

; NEIGHBORSEARCHING PARAMETERS
; nblist update frequency
nstlist                  = 5
; ns algorithm (simple or grid)
ns-type                  = grid
; Periodic boundary conditions: xyz (default), no (vacuum)
; or full (infinite systems only)
pbc                      = xyz
; nblist cut-off
rlist                    = 0.8
; OPTIONS FOR ELECTROSTATICS AND VDW
; Method for doing electrostatics
coulombtype              = Reaction-field
rcoulomb-switch          = 0
rcoulomb                 = 1.4
; Relative dielectric constant for the medium and the reaction field
epsilon_r                = 1
epsilon_rf               = """

	diec = str(dielectric)

	content2 = """
; Method for doing Van der Waals
vdw-type                 = cut-off
; cut-off lengths
rvdw-switch              = 0
rvdw                     = 1.4
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = no
; Extension of the potential lookup tables beyond the cut-off
table-extension          = 1
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; FFT grid size, when a value is 0 fourierspacing will be used
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters
pme_order                = 4
ewald_rtol               = 1e-05
ewald_geometry           = 3d
epsilon_surface          = 0
optimize_fft             = no


; Temperature coupling
Tcoupl          = berendsen
tc-grps         = System
tau_t           = 0.2
ref_t           = """

	temp = str(T)

	content3 = """

; Pressure coupling
Pcoupl                   = berendsen
Pcoupltype               = Isotropic
tau_p                    = 0.5
compressibility          = """

	compress = str(comp)

	content4 = """
ref_p                    = 1.0

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                 = no    ; Assign velocities to particles by taking them randomly from a Maxwell distribution

; OPTIONS FOR BONDS
constraints             = all-bonds
continuation            = yes
constraint_algorithm    = lincs
lincs_iter              = 1
lincs_order             = 4
"""

	inmd_h.write(content)
	inmd_h.write(diec)
	inmd_h.write(content2)
	inmd_h.write(temp)
	inmd_h.write(content3)
	inmd_h.write(compress)
	inmd_h.write(content4)
	inmd_h.close()

def main(ver,prefix,threshold):

	if not os.path.exists(prefix + ".LIQ.done"):
		if not os.path.exists(prefix + ".LIQ.cpt"):

			create_mdp(ver,prefix)

			gromppcommand = "grompp_" + ver + " -f in.LIQ.mdp -c " + prefix + ".NPT3.gro -p " + prefix + ".top -t " + prefix + ".NPT3.cpt -o " + prefix + ".LIQ.tpr -po " + prefix + ".LIQ.mdout.mdp -maxwarn 1"
			os.system(gromppcommand)

			os.system("mdrun_" + ver + " -deffnm " + prefix + ".LIQ")
			os.system('cp ' +prefix+'.LIQ.gro '+prefix+'.LIQ.done')

		else:
			os.system("mdrun_" + ver + " -deffnm " + prefix + ".LIQ -cpi " + prefix + ".LIQ.cpt -append")

print "\n>>> LIQ assemble alredy finished. Moving on...\n"

main(ver,prefix,threshold)
