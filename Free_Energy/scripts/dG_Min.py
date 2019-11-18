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


def Minimization1(lamb, moleculetype):
    inmin1 = "in.min1.mdp"
    inmin1_h = open(inmin1, 'w')
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
    ; No velocities during EM
    gen_vel                  = no
    ; options for bonds
    constraints              = all-bonds  ; we only have C-H bonds here
    ; Type of constraint algorithm
    constraint-algorithm     = lincs
    ; Do not constrain the starting configuration
    continuation             = no
    ; Highest order in the expansion of the constraint coupling matrix
    lincs-order              = 12

        """
    inmin1_h.write(content)
    inmin1_h.write(L)
    inmin1_h.write(content2)
    inmin1_h.write(mol)
    inmin1_h.write(content3)
    inmin1_h.close()

    if not os.path.exists(prefix + ".min1.gro"):
    	min1gromppcommand = "gmx_" + ver + " grompp -f in.min1.mdp -c " + prefix + ".min0.gro -p " + prefix + ".top -o " + prefix + ".min1.tpr -po " + prefix + ".min1.mdout.mdp -maxwarn 1"
    	os.system(min1gromppcommand)

    	min1command = "gmx_" + ver + " mdrun -deffnm " + prefix + ".min1"
    	os.system(min1command)
    else:
    	print "\n>>> Minimization 1 completed!\n"

def Minimization2(lamb, moleculetype):
    inmin2 = "in.min2.mdp"
    inmin2_h = open(inmin2, 'w')
    content = """; Run control
    integrator               = l-bfgs
    nsteps                   = 1000
    define                   = -DFLEXIBLE
    ; EM criteria and other stuff
    emtol                    = 10
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
    ; No velocities during EM
    gen_vel                  = no
    ; options for bonds
    constraints              = none     ; L-BFGS doesn't work with constraints
    ; Type of constraint algorithm
    constraint-algorithm     = lincs
    ; Do not constrain the starting configuration
    continuation             = no
    ; Highest order in the expansion of the constraint coupling matrix
    lincs-order              = 12

        """
    inmin2_h.write(content)
    inmin2_h.write(L)
    inmin2_h.write(content2)
    inmin2_h.write(mol)
    inmin2_h.write(content3)
    inmin2_h.close()

    if not os.path.exists(prefix + ".min2.gro"):
    	min1gromppcommand = "gmx_" + ver + " grompp -f in.min2.mdp -c " + prefix + ".min1.gro -p " + prefix + ".top -o " + prefix + ".min2.tpr -po " + prefix + ".min2.mdout.mdp -maxwarn 1"
    	os.system(min1gromppcommand)

    	min1command = "gmx_" + ver + " mdrun -deffnm " + prefix + ".min2"
    	os.system(min1command)
    else:
    	print "\n>>> Minimization 2 completed!\n"

########################################
Minimization1(lamb, molecule)
Minimization2(lamb, molecule)
