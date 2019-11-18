########################################
# Usage: 'python Equ.py ver prefix'
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

def create_mdp1(ver,prefix):
	############################################
	# equilibration NPT

	inNPT = "in.NPT1.mdp"
	inNPT_h = open(inNPT, 'w')
	content = """; VARIOUS PREPROCESSING OPTIONS
define          = -DFLEXIBLE

; RUN CONTROL PARAMETERS
integrator      = md
dt              = 0.002
nsteps          = 1500000

; OUTPUT CONTROL OPTIONS
nstxout         = 0
nstvout         = 0
nstenergy       = 100
nstlog          = 100
nstxtcout       = 0
xtc-precision   = 1000
energygrps      = System

; NEIGHBORSEARCHING PARAMETERS
nstlist         = 5
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
epsilon-r       = 1
epsilon-rf      = """

	diec = str(dielectric)

	content2 = """
Tcoupl          = berendsen
tc-grps         = System
tau_t           = 0.1
ref_t           = """
	temp = str(T)

	content3 = """
; Dispersion correction
DispCorr        = EnerPres

; Pressure coupling
Pcoupl          = berendsen
Pcoupltype      = Isotropic
tau_p           = 5.0
compressibility = """
	compress = str(comp)
	content4 = """
ref_p           = 100.0
refcoord_scaling = com

; GENERATE VELOCITIES FOR STARTUP RUN
gen_vel                 = yes    ; Assign velocities to particles by taking them randomly from a Maxwell distribution
gen_temp                = """

	content5 = """		; Temperature to generate corresponding Maxwell distribution
gen_seed                = -1     ; Seed for (semi) random number generation. Different numbers give different sets of velocities


; OPTIONS FOR BONDS
constraints     = all-bonds
continuation    = no
constraint_algorithm    = lincs
lincs_iter      = 1
lincs_order     = 4
		"""
	inNPT_h.write(content)
	inNPT_h.write(diec)
	inNPT_h.write(content2)
	inNPT_h.write(temp)
	inNPT_h.write(content3)
	inNPT_h.write(compress)
	inNPT_h.write(content4)
	inNPT_h.write(temp)
	inNPT_h.write(content5)
	inNPT_h.close()

def create_mdp2(ver,prefix):
	############################################
	# equilibration NPT

	inNPT = "in.NPT2.mdp"
	inNPT_h = open(inNPT, 'w')
	content = """; RUN CONTROL PARAMETERS
integrator               = md
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 1000000
; For exact run continuation or redoing part of a run
init_step                = 0
; mode for center of mass motion removal
comm-mode                = Linear
; number of steps for center of mass motion removal
nstcomm                  = 1
; group(s) for center of mass motion removal
comm-grps                =

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
; Output frequency for energies to log file and energy file
nstlog                   = 100
nstenergy                = 100
; Output frequency and precision for xtc file
xtc-precision            = 1000
; This selects the subset of atoms for the xtc file. You can
; select multiple groups. By default all atoms will be written.
xtc-grps                 =
; Selection of energy groups
energygrps               =

; NEIGHBORSEARCHING PARAMETERS
; nblist update frequency
nstlist                  = 5
; ns algorithm (simple or grid)
ns-type                  = grid
; Periodic boundary conditions: xyz (default), no (vacuum)
; or full (infinite systems only)
pbc                      = xyz
; nblist cut-off
rlist                    = 1.1

; OPTIONS FOR ELECTROSTATICS AND VDW
; Method for doing electrostatics
coulombtype              = PME
rcoulomb-switch          = 0
rcoulomb                 = 1.1
; Relative dielectric constant for the medium and the reaction field
epsilon_r                = 1
epsilon_rf               = """

	diec = str(dielectric)

	content2 = """
; Method for doing Van der Waals
vdw-type                 = cut-off
; cut-off lengths
rvdw-switch              = 0
rvdw                     = 1.1
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = no
; Extension of the potential lookup tables beyond the cut-off
table-extension          = 1
; Seperate tables between energy group pairs
energygrp_table          =
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

; GENERALIZED BORN ELECTROSTATICS
; Algorithm for calculating Born radii
gb_algorithm             = Still
; Frequency of calculating the Born radii inside rlist
nstgbradii               = 1
; Cutoff for Born radii calculation; the contribution from atoms
; between rlist and rgbradii is updated every nstlist steps
rgbradii                 = 2
; Salt concentration in M for Generalized Born models
gb_saltconc              = 0

; IMPLICIT SOLVENT (for use with Generalized Born electrostatics)
implicit_solvent         = No

; TEMPERATURE COUPLING
Tcoupl          = berendsen
tc-grps         = System
tau_t           = 0.1
ref_t           = """
	temp = str(T)

	content3 = """

; Pressure coupling
Pcoupl          = berendsen
Pcoupltype      = Isotropic
tau_p           = 5.0
compressibility = """
	compress = str(comp)
	content4 = """
ref_p           = 100.0
refcoord_scaling = com


; OPTIONS FOR BONDS
constraints     = all-bonds
continuation    = yes
constraint_algorithm    = lincs
lincs_iter      = 1
lincs_order     = 4
		"""
	inNPT_h.write(content)
	inNPT_h.write(diec)
	inNPT_h.write(content2)
	inNPT_h.write(temp)
	inNPT_h.write(content3)
	inNPT_h.write(compress)
	inNPT_h.write(content4)
	inNPT_h.close()

def create_mdp3(ver,prefix):
	############################################
	# equilibration NPT

	inNPT = "in.NPT3.mdp"
	inNPT_h = open(inNPT, 'w')
	content = """; RUN CONTROL PARAMETERS
integrator               = md
; Start time and timestep in ps
tinit                    = 0
dt                       = 0.002
nsteps                   = 1000000
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
nstxtcout                = 500
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

; GENERALIZED BORN ELECTROSTATICS
; Algorithm for calculating Born radii
gb_algorithm             = Still
; Frequency of calculating the Born radii inside rlist
nstgbradii               = 1
; Cutoff for Born radii calculation; the contribution from atoms
; between rlist and rgbradii is updated every nstlist steps
rgbradii                 = 2
; Salt concentration in M for Generalized Born models
gb_saltconc              = 0

; IMPLICIT SOLVENT (for use with Generalized Born electrostatics)
implicit_solvent         = No

; TEMPERATURE COUPLING
Tcoupl          = berendsen
tc-grps         = System
tau_t           = 0.1
ref_t           = """
	temp = str(T)

	content3 = """

; Pressure coupling
Pcoupl          = berendsen
Pcoupltype      = Isotropic
tau_p           = 0.5
compressibility = """
	compress = str(comp)
	content4 = """
ref_p           = 1.0
refcoord_scaling = com


; OPTIONS FOR BONDS
constraints     = all-bonds
continuation    = yes
constraint_algorithm    = lincs
lincs_iter      = 1
lincs_order     = 4
		"""
	inNPT_h.write(content)
	inNPT_h.write(diec)
	inNPT_h.write(content2)
	inNPT_h.write(temp)
	inNPT_h.write(content3)
	inNPT_h.write(compress)
	inNPT_h.write(content4)
	inNPT_h.close()

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
				   if line.split()[1] == "125":
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

	os.system("rm " + prefix + ".top")
 	os.system("mv " + prefix + ".out.top " + prefix + ".top")

def scaling(ver,prefix):
	########################################
	# Scaling system to 1000 molecules
	scalingcommand = "genconf_" + ver + " -f " + prefix + ".NPT2.gro -o " + prefix + ".1000.gro -nbox 2 2 2"
 	os.system(scalingcommand)

def drift_read(ver,prefix):
	####################################
	# Analyzing value of the Total-Energy drift
	os.system("echo Total-energy | g_energy_" + ver + " -s " + prefix + ".NPT3.tpr -f " + prefix + ".NPT3.edr -nmol " + str(N) + " -o " + prefix + ".equ-total_energy.xvg > " + prefix + ".equ-drift.out")

	f = open(prefix + ".equ-drift.out","r")
	string = "text"
	while True:
		string = f.readline()
		if not string: break
		if len(string.split())>0:
			if string.split()[0] == "Total":
				Energy = float(string.split()[5])
				break

	f2 = open(prefix + ".equ-drift.out","r")
	time = "text"
	while True:
		time = f2.readline()
		if not time: break
		if len(time.split())>0:
			if time.split()[0] == "Statistics":
				Time = float(time.split()[7])
				break
	return Energy,Time

def drift_output(ver,prefix,DF,name,T):
	ET = drift_read(ver,prefix)

	drift = abs(ET[0]*1000) # J/mol
	stime = ((ET[1])/1000) # ns

	value = drift/(DF*stime)

	print "####################################################################\n"
	print "Set of Parameters:"
	print "Simulation Time =", str(stime),"ns"
	print "Molecule name =", str(name)
	print "Degree of Freedom =", str(DF)
	print "Temperature =", str(T)

	return value


if not os.path.exists(prefix + ".NPT1.gro"):
   if not os.path.exists(prefix + ".NPT1.cpt"):

      create_mdp1(ver,prefix)

      gromppcommand = "grompp_" + ver + " -f in.NPT1.mdp -c " + prefix + ".min1.gro -p " + prefix + ".top -o " + prefix + ".NPT1.tpr -po " + prefix + ".NPT1.mdout.mdp -maxwarn 1"
      os.system(gromppcommand)

      mdruncommand = "mdrun_" + ver + " -nt 2 -deffnm " + prefix + ".NPT1"
      os.system(mdruncommand)
   else:
      mdruncommand = "mdrun_" + ver + " -nt 2 -deffnm " + prefix + ".NPT1 -cpi " + prefix + ".NPT1.cpt -append"
      os.system(mdruncommand)
else:
   print "\n>>> NPT1 assemble alredy finished. Moving on...\n"

if not os.path.exists(prefix + ".NPT2.gro"):
   if not os.path.exists(prefix + ".NPT2.cpt"):

      create_mdp2(ver,prefix)

      gromppcommand = "grompp_" + ver + " -f in.NPT2.mdp -c " + prefix + ".NPT1.gro -p " + prefix + ".top -t " + prefix + ".NPT1.cpt -o " + prefix + ".NPT2.tpr -po " + prefix + ".NPT2.mdout.mdp -maxwarn 1"
      os.system(gromppcommand)

      mdruncommand = "mdrun_" + ver + " -deffnm " + prefix + ".NPT2"
      os.system(mdruncommand)
   else:
      mdruncommand = "mdrun_" + ver + " -deffnm " + prefix + ".NPT2 -cpi " + prefix + ".NPT2.cpt -append"
      os.system(mdruncommand)
else:
   print "\n>>> NPT2 assemble alredy finished. Moving on...\n"

flag = 1

if not os.path.exists(prefix + ".NPT3.done"):
	out = "conv-eq.dat"
	out_t = open(out, 'w')
	content = 'lala'
	out_t.close()

	if not os.path.exists(prefix + ".NPT3.cpt"):
		if not os.path.exists(prefix + ".1000.gro"):
			scaling(ver,prefix)

		change_line(prefix + ".top", 1000)

		create_mdp3(ver,prefix)

		gromppcommand = "grompp_" + ver + " -f in.NPT3.mdp -c " + prefix + ".1000.gro -p " + prefix + ".top -o " + prefix + ".NPT3.tpr -po " + prefix + ".NPT3.mdout.mdp -maxwarn 1"
		os.system(gromppcommand)

		mdruncommand = "mdrun_" + ver + " -deffnm " + prefix + ".NPT3"
		os.system(mdruncommand)


		while flag == 1:
			out_t = open(out, 'a')
			value = drift_output(ver,prefix,DF,name,T)
			if abs(float(value))> threshold:

			  content = "\n >>> The calculated total drift of " + str(value) + " J/mol/ns/DegFree is higher than the previous stablished value of " + str(threshold) + " J/mol/ns/DegFree."
			  out_t.write(content)
			  out_t.close()
			  flag = 1

			  os.system("gmx_" + ver + " convert-tpr -s " + prefix + ".NPT3.tpr -extend 1000 -o " + prefix + ".NPT3.tpr")
			  os.system("rm *NPT3.tpr.* *NPT3.gro.* *drift.out *xvg")
			  os.system("mdrun_" + ver + " -deffnm " + prefix + ".NPT3 -cpi " + prefix + ".NPT3.cpt -append")

			else:
			  content = "\n >>> The calculated total drift is " + str(value) + " J/mol/ns/DegFree. This means that your simulation has converged (Threshold of " +str(threshold) +" J/mol/ns/DegFree). Congratulations! =)"
			  out_t.write(content)
			  out_t.close()
			  flag = 0
			  os.system('cp ' +prefix+'.NPT3.gro '+prefix+'.NPT3.done')
	else:
		mdruncommand = "mdrun_" + ver + " -deffnm " + prefix + ".NPT3 -cpi " + prefix + ".NPT3.cpt -append"
		os.system(mdruncommand)

		while flag == 1:
			out_t = open(out, 'a')
			value = drift_output(ver,prefix,DF,name,T)
			if abs(float(value))> threshold:

			  content = "\n >>> The calculated total drift of " + str(value) + " J/mol/ns/DegFree is higher than the previous stablished value of " + str(threshold) +" J/mol/ns/DegFree."
			  out_t.write(content)
			  out_t.close()
			  flag = 1

			  os.system("gmx_" + ver + " convert-tpr -s " + prefix + ".NPT3.tpr -extend 1000 -o " + prefix + ".NPT3.tpr")
			  os.system("rm *NPT3.tpr.* *NPT3.gro.* *drift.out *xvg")
			  os.system("mdrun_" + ver + " -deffnm " + prefix + ".NPT3 -cpi " + prefix + ".NPT3.cpt -append")

			else:
			  content = "\n >>> The calculated total drift is " + str(value) + " J/mol/ns/DegFree. This means that your simulation has converged (Threshold of " +str(threshold) +" J/mol/ns/DegFree). Congratulations! =)"
			  out_t.write(content)
			  out_t.close()
			  flag = 0
			  os.system('cp ' +prefix+'.NPT3.gro '+prefix+'.NPT3.done')
else:
   print "\n>>> NPT3 assemble alredy finished. Moving on...\n"
