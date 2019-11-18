import os
import sys

ver = sys.argv[1]
prefix = sys.argv[2]
lamb = sys.argv[3]

def grep_molecule(file1):
	f2 = open(file1,"r")
	l2 = "text"

	while True:
		l2 = f2.readline()
		if not l2: break
		if len(l2.split()) > 0:
			if l2[2:10] == 'Compound':
				l2 = f2.readline()
				molecule = l2.split()[0]

	return molecule

molecule = grep_molecule(prefix+".top")

def Equilibration1(lamb,moleculetype):
    inequ1 = "in.NVT.mdp"
    inequ1_h = open(inequ1, 'w')
    content = """; Run control
    integrator               = sd       ; Langevin dynamics
    tinit                    = 0
    dt                       = 0.002
    nsteps                   = 50000    ; 100 ps
    nstcomm                  = 100
    ; Output control
    nstxout                  = 0
    nstvout                  = 0
    nstfout                  = 0
    nstlog                   = 500
    nstenergy                = 500
    nstxout-compressed       = 500
  
    ; Neighborsearching and short-range nonbonded interactions
    cutoff-scheme            = verlet
    nstlist                  = 10
    ns_type                  = Grid
    pbc                      = xyz
    rlist                    = 0.8
    ; Electrostatics
    coulombtype              = Reaction-field
    rcoulomb                 = 1.4

    ; van der Waals
    vdwtype                  = cutoff
    rvdw                     = 1.4
    ; Apply long range dispersion corrections for Energy and Pressure
    DispCorr                  = EnerPres
    ; Spacing for the PME/PPPM FFT grid
    fourierspacing           = 0.12
    ; EWALD/PME/PPPM parameters
    pme_order                = 6
    ewald_rtol               = 1e-06
    epsilon_surface          = 0
    ; Temperature coupling
    ; tcoupl is implicitly handled by the sd integrator
    tc_grps                  = system
    tau_t                    = 1.0
    ref_t                    = 300
    ; Pressure coupling is off for NVT
    Pcoupl                   = No
    tau_p                    = 0.5
    compressibility          = 4.5e-05
    ref_p                    = 1.0
    ;################################################################################
    ; Free energy control stuff
    free_energy              = yes
    init_lambda_state        = """

    L = str(lamb)

    content2 = """
    delta_lambda             = 0
    calc_lambda_neighbors    = 1        ; only immediate neighboring windows
    ; Vectors of lambda specified here
    ; Each combination is an index that is retrieved from init_lambda_state for each simulation
    ; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48   49
    vdw_lambdas              = 0.00 0.02 0.04 0.07 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.93 0.96 0.98 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    coul_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.02 0.04 0.07 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.93 0.96 0.98 1.00
    ; We are not transforming any bonded or restrained interactions
    bonded_lambdas           = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
    restraint_lambdas        = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
    ; Masses are not changing (particle identities are the same at lambda = 0 and lambda = 1)
    mass_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
    ; Not doing simulated temperting here
    temperature_lambdas      = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
    ; Options for the decoupling
    sc-alpha                 = 0.5
    sc-coul                  = no       ; linear interpolation of Coulomb (none in this case)
    sc-power                 = 1
    sc-sigma                 = 0.3
    couple-moltype           = """

    mol = str(molecule)

    content3 = """
    couple-lambda0           = none     ; only van der Waals interactions
    couple-lambda1           = vdw-q    ; turn off everything, in this case only vdW
    couple-intramol          = no
    nstdhdl                  = 10
    ;################################################################################
    ; Generate velocities to start
    gen_vel                  = yes
    gen_temp                 = 300
    gen_seed                 = 7
    ; options for bonds
    constraints              = all-bonds  ; we only have C-H bonds here
    ; Type of constraint algorithm
    constraint-algorithm     = lincs
    ; Do not constrain the starting configuration
    continuation             = no
    ; Highest order in the expansion of the constraint coupling matrix
    lincs-order              = 12
        """
    inequ1_h.write(content)
    inequ1_h.write(L)
    inequ1_h.write(content2)
    inequ1_h.write(mol)
    inequ1_h.write(content3)
    inequ1_h.close()

    if not os.path.exists(prefix + ".NVT.gro"):
    	min1gromppcommand = "gmx_" + ver + " grompp -f in.NVT.mdp -c " + prefix + ".min2.gro -p " + prefix + ".top -o " + prefix + ".NVT.tpr -po " + prefix + ".NVT.mdout.mdp -maxwarn 1"
    	os.system(min1gromppcommand)

    	min1command = "gmx_" + ver + " mdrun -deffnm " + prefix + ".NVT"
    	os.system(min1command)
    else:
    	print "\n>>> NVT completed!\n"


def Equilibration2(lamb,moleculetype):
    inequ2 = "in.NPT.mdp"
    inequ2_h = open(inequ2, 'w')
    content = """; Run control
    integrator               = sd       ; Langevin dynamics
    tinit                    = 0
    dt                       = 0.002
    nsteps                   = 50000  ; 100 ps
    nstcomm                  = 100
    ; Output control
    nstxout                  = 0
    nstvout                  = 0
    nstfout                  = 0
    nstlog                   = 500
    nstenergy                = 500
    nstxout-compressed       = 500
    ; Neighborsearching and short-range nonbonded interactions
    nstlist                  = 10
    ns-type                  = Grid
    pbc                      = xyz
    rlist                    = 0.8

    ; OPTIONS FOR ELECTROSTATICS AND VDW
    coulombtype              = Reaction-field
    pme_order                = 4    ; cubic interpolation
    fourierspacing           = 0.12 ; grid spacing for FFT
    rcoulomb                 = 1.4
    vdw-type                 = cut-off
    rvdw                     = 1.4

    ; Apply long range dispersion corrections for Energy and Pressure
    DispCorr                  = EnerPres
    ; EWALD/PME/PPPM parameters
    ewald_rtol               = 1e-06
    epsilon_surface          = 0
    ; Temperature coupling
    ; tcoupl is implicitly handled by the sd integrator
    tc_grps                  = system
    tau_t                    = 1.0
    ref_t                    = 300
    ; Pressure coupling is off for NVT
    Pcoupl                   = Parrinello-Rahman
    tau_p                    = 1.0
    compressibility          = 4.5e-05
    ref_p                    = 1.0
    ;################################################################################
    ; Free energy control stuff
    free_energy              = yes
    init_lambda_state        = """

    L = str(lamb)

    content2 = """
    delta_lambda             = 0
    calc_lambda_neighbors    = 1        ; only immediate neighboring windows
    ; Vectors of lambda specified here
    ; Each combination is an index that is retrieved from init_lambda_state for each simulation
    ; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48   49
    vdw_lambdas              = 0.00 0.02 0.04 0.07 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.93 0.96 0.98 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
    coul_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.02 0.04 0.07 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.93 0.96 0.98 1.00
    ; We are not transforming any bonded or restrained interactions
    bonded_lambdas           = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
    restraint_lambdas        = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
    ; Masses are not changing (particle identities are the same at lambda = 0 and lambda = 1)
    mass_lambdas             = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
    ; Not doing simulated temperting here
    temperature_lambdas      = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
    ; Options for the decoupling
    sc-alpha                 = 0.5
    sc-coul                  = no       ; linear interpolation of Coulomb (none in this case)
    sc-power                 = 1
    sc-sigma                 = 0.3
    couple-moltype           = """

    mol = str(molecule)

    content3 = """
    couple-lambda0           = none     ; only van der Waals interactions
    couple-lambda1           = vdw-q    ; turn off everything, in this case only vdW
    couple-intramol          = no
    nstdhdl                  = 10
    ;################################################################################
    ; Generate velocities to start
    gen_vel                  = no
    gen_temp                 = 300
    gen_seed                 = 7
    ; options for bonds
    constraints              = all-bonds  ; we only have C-H bonds here
    ; Type of constraint algorithm
    constraint-algorithm     = lincs
    ; Do not constrain the starting configuration
    continuation             = no
    ; Highest order in the expansion of the constraint coupling matrix
    lincs-order              = 12
        """
    inequ2_h.write(content)
    inequ2_h.write(L)
    inequ2_h.write(content2)
    inequ2_h.write(mol)
    inequ2_h.write(content3)
    inequ2_h.close()

    if not os.path.exists(prefix + ".NPT.gro"):
    	min1gromppcommand = "gmx_" + ver + " grompp -f in.NPT.mdp -c " + prefix + ".NVT.gro -p " + prefix + ".top -t " + prefix + ".NVT.cpt -o " + prefix + ".NPT.tpr -po " + prefix + ".NPT.mdout.mdp -maxwarn 1"
    	os.system(min1gromppcommand)

    	min1command = "gmx_" + ver + " mdrun -deffnm " + prefix + ".NPT"
    	os.system(min1command)
    else:
    	print "\n>>> NPT completed!\n"


Equilibration1(lamb,molecule)
Equilibration2(lamb,molecule)
