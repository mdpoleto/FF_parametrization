##################################################################################################################################
# Usage: 'python MM-Calc.py OFF prefix D1 D2 D3 D4'
# Example: 'python MM-Calc.py OFF ben C1 C2 C3 C4'
#
# This script was written to scan the energy for dihedrals every 30 degrees (0-360), following the steps:
#   1-  It assumes that you already have the needed PDB structures named 'prefix_dihedral' (e.g. ben_270.pdb, ben_300.pdb).
#   2-  It assumes that you have a force field (example "gromos53a6.ff" or any other name) in your working directory that will be required at 'pdb2gmx' step
#   3-  It will create a topology based on .rtp database. Make sure your topology is working there!
#   4-  If you choose "OFF", it will suppress the wanted 'gd_' to calculate the potential without the dihedral contribution.
#   5-  If you choose "ON", it will not suprees the wanted 'gd_', so it will calculate the total potential force.
#   6-  It will insert the 'dihedral_restraint' lines into the .top file for your each dihedral angle evaluated.
#   7-  It will run a BLFGS minimization with double precision. Make sure you have a working version of GROMACS! You can change it at line #30
#   8-  If you choose "OFF", it will grep the energy and write it on a prefix.MM.dat file.
#   9-  If you choose "ON", it will grep the energy and write it on a prefix.FIT.dat file.
#   10- All maximum forces will be printed on the convergence.dat file.
#
# (In case of fail, contact 'marcelodepolo@gmail.com'#############################################################################

#########################################
# Loading libraries
import os
import sys
import math
import operator
import linecache
from pymol import cmd

#######################################
# Loading prefix arguments
ver = '507_d'       # gromacs version (e.g. 504)
prefix = sys.argv[2]    # molecule name (e.g. ben)

#######################################
# Loading Dihedral Atoms names

D1 = str(sys.argv[3])
D2 = str(sys.argv[4])
D3 = str(sys.argv[5])
D4 = str(sys.argv[6])

######################################

fit_flag = sys.argv[1] #'ON' # or OFF



def create_mdp(ver):
	############################################
	# Conjugate Gradient energy minimization

	inmin = "in.min.mdp"
	inmin_h = open(inmin, 'w')
	content = """; Define can be used to control processes

			title		    =  Yo
			cpp                 =  /lib/cpp
			constraints         =  none
			integrator          =  l-bfgs
			emtol		    =  0.01
			emstep              =  0.1
			nstcgsteep	    =  100
			nsteps		    =  500000
			nstcomm             =  1
			nstxout             =  100
			nstvout             =  100
			nstfout             =  0
			nstlog              =  100
			nstenergy           =  100
			nstlist             =  10
			ns_type             =  grid
			coulombtype         =  PME
			vdw-type			=  PME
			rlist               =  0.9
			rcoulomb            =  0.9
			rvdw                =  0.9
			pme_order           =  4
			; Berendsen temperature coupling is on in four groups
			Tcoupl              =  no
			; Isotropic pressure coupling is now on
			Pcoupl              =  no
			dihre = yes
			dihre_fc = 100000
			dihre_tau = 0
			nstdihreout = 1 """

	inmin_h.write(content)
	inmin_h.close()

def create_topol_box(ver,name):

	###########################################
	# setup periodic box
	pdb2gmxcommand = "gmx_" + ver + " pdb2gmx -ff gromos53a6 -water spc -f " + name + ".pdb -p " + name + ".top -o " + name + ".gro -i " + name + ".posre.itp"
	os.system(pdb2gmxcommand)

	###########################################
	# setup periodic box
	setupcommand = "gmx_" + ver + " editconf -f " + name + ".gro -o " + name + ".pbc.gro -c -d 1.0"
	os.system(setupcommand)

def atoms_label(name):

		f = open(name + ".top","r")
		l = ""
		dic = {}

		while l[:9]!= "; residue":		#acessing the atoms
			l = f.readline()
			if not l: break

		while l[:9]!= "[ bonds ]":		#extraction the atoms
			l = f.readline()
			if not l: break
			if  len(l.split())==11:
				atom_number =  l.split()[0]
				atom_name =  l.split()[4]
				dic[atom_name] = atom_number
		f.close()

		return dic

def change_line(file,dic):

		f_in = open(file,"r")
		f_out = open(file[:-4] + ".out.top","w")
		line = "anygiventext"
		new_line = ""
		while True:
			line = f_in.readline()
			if not line: break
			if len(line.split()) > 0:
				   if line.split()[0] == dic[D1]:
					   if line.split()[1] == dic[D2]:
						   if line.split()[2] == dic[D3]:
							   if line.split()[3] == dic[D4]:
								   new_line = ";" + line
								   f_out.write(new_line)
							   else:
								   f_out.write(line)
						   else:
							   f_out.write(line)
					   else:
						   f_out.write(line)

				   elif line.split()[0] == dic[D4]:
					   if line.split()[1] == dic[D3]:
						   if line.split()[2] == dic[D2]:
							   if line.split()[3] == dic[D1]:
								   new_line = ";" + line
								   f_out.write(new_line)
							   else:
								   f_out.write(line)
						   else:
							   f_out.write(line)
					   else:
						   f_out.write(line)
				   else:
					f_out.write(line)
			else:
				f_out.write(line)
		f_in.close()
		f_out.close()
		os.system("rm " + file)
		os.system("mv " + file[:-4]+ ".out.top " + file)

def insert_line(file,dic,i):

		f_in = open(file,"r")
		f_out = open(file[:-4] + ".out.top","w")
		line = "anygiventext"
		new_line = ""

		while True:
			line = f_in.readline()
			if not line: break
			if len(line.split()) > 0:
				   if line.split()[0] == "[":
					   if line.split()[1] == "angles":
						   new_line = "[ dihedral_restraints ] \n; ai   aj    ak    al  type   phi   dphi  kfac\n"
						   atoms_D = "  " + dic[D1] + "   " + dic[D2] + "    " + dic[D3] + "    " + dic[D4] + "    1   " + str(i) + "   0    10000\n"
						   angle = "\n[ angles ]\n"
						   f_out.write(new_line)
						   f_out.write(atoms_D)
						   f_out.write(angle)
					   else:
						   f_out.write(line)
				   else:
					f_out.write(line)
			else:
				f_out.write(line)
		f_in.close()
		f_out.close()

		os.system("rm " + file)
		os.system("mv " + file[:-4]+ ".out.top " + file)

def calc_ener(ver,name):
	###########################################
	# Calculating Energies

	gromppcommand = "gmx_" + ver + " grompp -f in.min.mdp -c " + name + ".pbc.gro -p " + name + ".top -o " + name + ".MIN.tpr -po " + name + ".MIN.mdout.mdp" " -maxwarn 1"
	os.system(gromppcommand)

	mdruncommand = "gmx_" + ver + " mdrun -deffnm " + name + ".MIN" " -nt 2 -c " + name + ".MIN.pdb"
	os.system(mdruncommand)

	os.system("echo Potential | gmx_" + ver + " energy -f " + name + ".MIN.edr -s " + name + ".MIN.tpr -o " + name + ".pot.xvg")

def grep_value(prefix,name,i,dic_energy):

	def file_len(fname):
		with open(fname) as f:
			for i, l in enumerate(f):
				pass
		return i + 1

	num_lines = file_len(name + ".pot.xvg")
	linha = linecache.getline(name + ".pot.xvg", num_lines)
	potential =  linha.split()[1]

	energy = float(potential)

	dic_energy[i] = energy

	return energy


def grep_force(name):

	f = open(name + ".MIN.log", 'r')
	line = f.readline()

	while True:
		line = f.readline()
		if not line: break
		if len(line.split()) > 0:
			if line.split()[0] == 'Maximum' and line.split()[1] == 'force':
				force = float(line.split()[3])

	return force

def main(ver,prefix,D1,D2,D3,D4):
	#########################################
	# Setting Degree step

	degrees = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360]
	dic ={}
	dic_energy ={}
	dic_normalized ={}

	conv = open('convergence.dat', 'w')
	conv.write("Degree      Max.Force(kJ/mol)\n")

	create_mdp(ver)

	for i in degrees:
		name = prefix + "_" + str(i)

		create_topol_box(ver,name)
		if i == 0:
			dic = atoms_label(name)

			os.system("cp " + name + ".top topology_check.top")

		if fit_flag == 'OFF':
			change_line(name + ".top",dic)
		else:
			pass

		insert_line(name + ".top", dic,i)

		calc_ener(ver,name)
		if i == 0:
			os.system("cp " + name + ".gro structure_min_check.gro")
		grep_value(prefix,name,i,dic_energy)


		force = grep_force(name)
		deg = str(i).ljust(14,' ')
		conv.write(deg)
		conv.write(str(force))
		conv.write('\n')

	conv.close()


	if fit_flag == 'OFF':
		f_out = open(prefix + ".MM.dat","w")
	else:
		f_out = open(prefix + ".FIT.dat","w")


	minimum = min(dic_energy.values())

	for k,v in dic_energy.items():
		dic_normalized[k] = v - minimum
	for k,v in sorted(dic_normalized.items()):
		string = str(k) + " " + str(v) + "\n"
		f_out.write(string)
	f_out.close()



	print "##################################################"
	print "\n\nWe have found the following atoms:\n"
	print dic ,"\n"


	print "\n>>> Removing useless files..."
	os.system("rm *itp " + prefix + "*.gro *trr *edr *tpr *log " + prefix + "*.top *xvg *_*.MIN.pdb *mdp")
	print ">>> Done!\n"

	print "##################################################"
	print "\n>>> Hey! I've finished to calculate MM Energies for your dihedral! Check the .dat output!"
	print "\n>>> We assumed that your input structure and topology (on .rtp database) were correct.\n If you have any doubts, please check structure_check.gro and topology_check.top."
	print "\n>>> See ya! =)\n"

main(ver,prefix,D1,D2,D3,D4)
