########################################
# Usage: 'python Analysis-PropLiq.py ver prefix'

# This script was written to calculate Liquids properties and compare their values with experimental ones. For that, it assumes some folder and file structures:
#
# I   - LIQ files with the format 'prefix.NPTx.[edr,tpr,xtc,etc...]'
# II  - SD files with the format 'prefix.SIM.[edr,tpr,xtc,etc...]'
#
# All simulation scripts can be found by the names Prep.py, Min.py, Equ.py, LIQ.py, SD.py and DOS.py. You can see the protocol by reading them.

#########################################
# Loading libraries
import os
import sys
import math
import numpy
from scipy import stats

#######################################
# Loading prefix arguments
ver = sys.argv[1]       # gromacs version (e.g. 507)
prefix = sys.argv[2]    # molecule name (e.g. ben)
ver2 = '506'
#######################################
# Constant values
Kb = 0.0083 # kJ/mol/K
N = 1000 # Number of molecules in system
threshold = 0.5 # J/mol/ns/DegFree (Total-Energy drift)
eA = 4.80320440079 # Debye

Cvib = 0

###############################################
# Reading Molecules Library (DF = Degree of freedom; ExD = Experimental Density in g/cm^3; ExV = Experimental Vaporization Enthalpy kJ/mol; ExC = Isothermal Compressibility in 1/GPa; ExE = Thermal Expansion Coefficient in 10^-3/K)

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

Tbelow = T-20
Tabove = T+20

DF = int(prop[3])
if prop[4] == 'nan':
	ExD = prop[4]
else:
	ExD = float(prop[4]) # g/cm^3

if prop[5] == 'nan':
	ExV = prop[5]
else:
	ExV = float(prop[5]) # kJ/mol

if prop[6] == 'nan':
	ExE = prop[6]
else:
	ExE = float(prop[6]) # 10^-3/K

if prop[7] == 'nan':
	ExC = prop[7]
else:
	ExC = float(prop[7]) # 1/GPa

if prop[8] == 'nan':
	ExDc = prop[8]
else:
	ExDc = float(prop[8])

if prop[9] == 'nan':
	ExDm = prop[9]
else:
	ExDm = float(prop[9]) # D

if prop[10] == 'nan':
	ExCp = prop[10]
else:
	ExCp = float(prop[10]) # J/mol*K

ExCv = ExCp

if prop[12] == 'nan':
	MW = prop[12]
else:
	MW = float(prop[12]) # g/mol

###################################################################
# Setting some variables

block = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000]

fout = open("TDprops.txt","w")
header = '##################################################################\n#Thermodynamical Properties of Liquid ' + str(name) + '\n'
fout.write(header)

def calc_dipole(file1,file2):

	f2 = open(file2,"r")
	l2 = "text"
	f1 = open(file1,"r")
	l1 = "text"
	nl = len(f1.readlines())
	f1.seek(0)
	dic1 = {}
	dic2 = {}
	dic3 = {}


	while l2[2:9] != "residue":		#acessing the atoms
		l2 = f2.readline()
		if not l2: break

	while l2[2:7] != 'bonds':		#extraction the atoms
		l2 = f2.readline()
		if not l2: break
		if len(l2.split()) > 1:
			if l2.split()[0] == '[': break
		if len(l2) > 1:
			atom_code = l2.split()[4]
			atom_charge = float(l2.split()[6])
			dic2[atom_code] = atom_charge

	for i in range(2):
		f1.readline()
	for i in range(nl-3):
		l1 = f1.readline()
		atom_code = l1.split()[1]
		atom_X = float(l1.split()[3])*10
		atom_Y = float(l1.split()[4])*10
		atom_Z = float(l1.split()[5])*10
		if atom_code in dic2.keys():
			dic1[atom_code] = (atom_X*float(dic2[atom_code]), atom_Y*float(dic2[atom_code]), atom_Z*float(dic2[atom_code]))
			dic3[atom_code] = (atom_X, atom_Y, atom_Z)

	f1.close()
	f2.close()

	listX = []
	listY = []
	listZ = []

	for i in dic1:

		compX = dic1[i][0]
		listX.append(compX)

		compY = dic1[i][1]
		listY.append(compY)

		compZ = dic1[i][2]
		listZ.append(compZ)


	dipoleX = (numpy.sum(listX))*eA
	dipoleY = (numpy.sum(listY))*eA
	dipoleZ = (numpy.sum(listZ))*eA

	module = round(((dipoleX**2) + (dipoleY**2) + (dipoleZ**2))**0.5,2)

	return module

def calc_prop(ver,prefix, block):
	###########################################
	# Creating folders and input

	os.system("mkdir dHvap")
	os.system("mkdir Density")
	os.system("mkdir TexpC")
	os.system("mkdir Comp")
	os.system("mkdir MSD")
	os.system("mkdir Cv")
	os.system("mkdir Diec")

	##############################################
	# Calculating dipole for topology used
	dipole = calc_dipole(prefix + '.gro',prefix + '.top')

	##############################################
	# Calculating Total-energy for drift
	os.system("echo Total-energy | g_energy_" + ver + " -s " + prefix + ".NPT3.tpr -f " + prefix + ".NPT3.edr -nmol " + str(N) + " -o " + prefix + ".equ-total_energy.xvg > " + prefix + ".equ-drift.out")
	os.system("echo Total-energy | g_energy_" + ver + " -s " + prefix + ".LIQ.tpr -f " + prefix + ".LIQ.edr -nmol " + str(N) + " -o " + prefix + ".LIQ-total_energy.xvg")

	f = open(prefix + ".LIQ-total_energy.xvg","r")
	string = "text"

	string = f.readlines()[-1]
	line = string.strip('\n')
	time_ps = int(round(float(line.split()[0])))

	##############################################
	# Ensemble average dependent properties

	# Density
	for t in block:
		os.system("echo Density | g_energy_" + ver + " -s " + prefix + ".LIQ.tpr -f " + prefix + ".LIQ.edr -b " + str(t) +" -e "+str(t+1000) + " -nmol " + str(N) + " -o Density/" + prefix + ".LIQ.density-"+str(t)+".xvg > Density/" + prefix + ".LIQ.density-"+str(t)+".out")

	# Liquid's diffusion
	os.system("echo System | g_msd_" + ver + " -s " + prefix + ".NPT3.gro -f " + prefix + ".NPT3.xtc -b 500 -e 1000 -o MSD/" + prefix + ".NPT3.beg-msd.xvg > MSD/" + prefix + ".beg-msd.out")
	os.system("echo System | g_msd_" + ver + " -s " + prefix + ".NPT3.gro -f " + prefix + ".NPT3.xtc -b 1500 -e 2000 -o MSD/" + prefix + ".NPT3.end-msd.xvg > MSD/" + prefix + ".end-msd.out")

	##############################################
	# Fluctations dependent properties / block averaging

	# Heat of vaporization
	for t in block:
		os.system("echo Potential | g_energy_" + ver + " -s " + prefix + ".LIQ.tpr -f " + prefix + ".LIQ.edr -b " + str(t) +" -e "+str(t+1000)+" -o dHvap/" + prefix + ".LIQ.potential-"+str(t)+".xvg > dHvap/" + prefix + ".epotsys-"+str(t)+".out")
		os.system("echo Potential | g_energy_" + ver + " -s " + prefix + ".SD.tpr -f " + prefix + ".SD.edr -b " + str(t*10) + " -e "+str((t*10)+10000)+" -o dHvap/" + prefix + ".SD.potential-"+str(t*10)+".xvg > dHvap/" + prefix + ".epotg-"+str(t*10)+".out")

	# Dielectric Constant
	dipole = calc_dipole(prefix + '.gro',prefix + '.top')

	#for i in range(0,time_ps,1000):
		#os.system("echo System | g_dipoles_" + ver + " -s " + prefix + ".LIQ.tpr -f " + prefix + ".LIQ.xtc -b 0 -e "+str(i+1000) + " -corr mol -mu " + str(ExDm) + " -temp " + str(T) + " -mumax 5 -o Diec/" + prefix + ".Mu-"+str(i)+".xvg -eps Diec/" + prefix + ".Diec-"+str(i)+".xvg -d Diec/" + prefix + ".dipdist-"+str(i)+".xvg -c Diec/" + prefix + ".dipcorr-"+str(i)+".xvg -a Diec/" + prefix + ".aver-"+str(i)+".xvg > Diec/" + prefix + ".Diec-"+str(i)+".out ")

	# Heat Capacity
	os.system("echo Total-energy | g_energy_" + ver + " -s " + prefix + ".T0.tpr -f " + prefix + ".T0.edr -b 1000  -nmol " + str(N) + " -o Cv/" + prefix + ".Etot-T0.xvg > Cv/" + prefix + ".Etot-T0.out")
	os.system("echo Total-energy | g_energy_" + ver + " -s " + prefix + ".Tlow.tpr -f " + prefix + ".Tlow.edr -b 1000 -nmol " + str(N) + " -o Cv/" + prefix + ".Etot-Tlow.xvg > Cv/" + prefix + ".Etot-Tlow.out")
	os.system("echo Total-energy | g_energy_" + ver + " -s " + prefix + ".Thigh.tpr -f " + prefix + ".Thigh.edr -b 1000 -nmol " + str(N) + " -o Cv/" + prefix + ".Etot-Thigh.xvg > Cv/" + prefix + ".Etot-Thigh.out")

	# Thermal Expansion Coefficient
	os.system("echo Density | g_energy_" + ver + " -s " + prefix + ".T0.tpr -f " + prefix + ".T0.edr -b 1000 -nmol " + str(N) + " -o TexpC/" + prefix + ".Den-T0.xvg > TexpC/" + prefix + ".Den-T0.out")
	os.system("echo Density | g_energy_" + ver + " -s " + prefix + ".Tlow.tpr -f " + prefix + ".Tlow.edr -b 1000 -nmol " + str(N) + " -o TexpC/" + prefix + ".Den-Tlow.xvg > TexpC/" + prefix + ".Den-Tlow.out")
	os.system("echo Density | g_energy_" + ver + " -s " + prefix + ".Thigh.tpr -f " + prefix + ".Thigh.edr -b 1000 -nmol " + str(N) + " -o TexpC/" + prefix + ".Den-Thigh.xvg > TexpC/" + prefix + ".Den-Thigh.out")

	# Isothermal Compressibility
	os.system("echo Pressure | g_energy_" + ver + " -s " + prefix + ".P0.tpr -f " + prefix + ".P0.edr -b 1000 -nmol " + str(N) + " -o Comp/" + prefix + ".Press-P0.xvg > Comp/" + prefix + ".Press-P0.out")
	os.system("echo Pressure | g_energy_" + ver + " -s " + prefix + ".Plow.tpr -f " + prefix + ".Plow.edr -b 1000 -nmol " + str(N) + " -o Comp/" + prefix + ".Press-Plow.xvg > Comp/" + prefix + ".Press-Plow.out")
	os.system("echo Pressure | g_energy_" + ver + " -s " + prefix + ".Phigh.tpr -f " + prefix + ".Phigh.edr -b 1000 -nmol " + str(N) + " -o Comp/" + prefix + ".Press-Phigh.xvg > Comp/" + prefix + ".Press-Phigh.out")

def drift_read(ver,prefix):
	####################################
	# Analyzing value of the Total-Energy drift
	t = open(prefix + ".LIQ-total_energy.xvg","r")
	string = "text"

	string = t.readlines()[-1]
	line = string.strip('\n')
	time_ps = int(round(float(line.split()[0])))


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
	return Energy,Time, time_ps

def drift_output(ver,prefix,DF,name,T):
	ET = drift_read(ver,prefix)

	drift = abs(ET[0]*1000) # J/mol
	stime = ((ET[1])/1000) # ns

	value = drift/(DF*stime) # J/mol/ns/DegFree

	print "####################################################################"
	print "Set of Parameters:"
	print "> Convergence Time needed =", str(stime),"ns"
	print "> Total drift energy =", str(round(value,2)), "J/mol/ns/DegFree"
	print "> Molecule name =", str(name)
	print "> Degree of Freedom =", str(DF)
	print "> Temperature =", str(T)

	return value

def calc_diffusion(ver, prefix):
	#####################################################
	# Colecting Initial Diffusion of Liquid system (Dif)
	f3 = open("MSD/" +prefix + ".beg-msd.out","r")
	dif1 = "text"
	while True:
		  dif1 = f3.readline()
		  if not dif1: break
		  if len(dif1.split())>0:
			 if dif1.split()[0] == "D[":
				break

	Difbeg = float(dif1.split()[2])

	#####################################################
	# Colecting Final Diffusion of Liquid system (Dif)
	f4 = open("MSD/" +prefix + ".end-msd.out","r")
	dif2 = "text"
	while True:
		  dif2 = f4.readline()
		  if not dif2: break
		  if len(dif2.split())>0:
			 if dif2.split()[0] == "D[":
				break

	Difend = float(dif2.split()[2])

	DeltaD = 2*(Difend - Difbeg)/(Difend + Difbeg)

	if abs(float(DeltaD)) > 0:
		print ">>> |DeltaD| => "+ str(round(abs(DeltaD),4))+ " > 0, it means your system did not freeze over the simulation.\n"
	elif abs(float(DeltaD)) > 0.5:
		print ">>> |DeltaD| => "+ str(round(abs(DeltaD),4))+ " > 0.5, it seams that your system froze over the simulation. Check it manually.\n"

	string0 = '\n|DeltaD| = ' + str(abs(DeltaD)) + '\n'
	fout.write(string0)

	return round(DeltaD,4)

def calc_den(ver,prefix):
	#####################################################
	# Colecting Density of Liquid system (Den)

	list_den = []

	for t in block:
		f0 = open("Density/"+ prefix + ".LIQ.density-"+str(t)+".out","r")
		den = "text"
		while True:
			den = f0.readline()
			if not den: break
			if len(den.split())>0:
				if den.split()[0] == "Density":
					break

		density = float(den.split()[1])/1000 #convert to g/cm^3
		list_den.append(density)
	Density = round(float(numpy.average(list_den)),4)
	Std_Dev = round(float(numpy.std(list_den)),2)

	if not ExD == 'nan':
		Error1 = round(abs((ExD - Density)/ExD)*100,4)
	else:
		Error1 = 'nan'

	print ">>> Aver. Density =", str(Density), "+-", str(Std_Dev), "g/cm^3"
	print ">>> Experimental =", ExD,"g/cm^3"
	print ">>> Error =", Error1,"%\n"

	string1 = '\nDensity = ' + str(Density) + ' +- ' + str(Std_Dev) + ' g/cm^3'
	string1_1 = '\nExperimental = ' + str(ExD) + ' g/cm^3'
	string1_2 = '\nError = ' + str(Error1) + '%\n'
	fout.write(string1)
	fout.write(string1_1)
	fout.write(string1_2)

	return Density, Error1

def calc_dHvap(ver,prefix,block):
	#####################################################
	# Colecting Potential energy of Gas system (Epotg)
	list_Epotg = []
	for t in block:
		f1 = open("dHvap/" + prefix + ".epotg-"+str(t*10)+".out","r")
		string1 = "text"
		while True:
			string1 = f1.readline()
			if not string1: break
			if len(string1.split())>0:
				if string1.split()[0] == "Potential":
					break
		Epotg = float(string1.split()[1])
		list_Epotg.append(Epotg)

	#####################################################
	# Colecting Potential energy of Liquid system (Epotl)
	list_Epotl = []
	for t in block:
		f2 = open("dHvap/" + prefix + ".epotsys-"+str(t)+".out","r")
		string2 = "text"
		while True:
			string2 = f2.readline()
			if not string2: break
			if len(string2.split())>0:
				if string2.split()[0] == "Potential":
					break

		Epotl = float(string2.split()[1])/N
		list_Epotl.append(Epotl)

	list_dHvap = []

	for i in range(len(list_Epotl)):
		Epotg = float(list_Epotg[i])
		Epotl = float(list_Epotl[i])
		dHvap = (Epotg + (Kb*298) - Epotl) # originally it is Kb*T, but the reference temperature is 298 K.
		list_dHvap.append(dHvap)

	Enthalpy = round(float(numpy.average(list_dHvap)),4)
	Std_Dev = round(float(numpy.std(list_dHvap)),2)

	if not ExV == 'nan':
		Error2 = round(abs((ExV - dHvap)/ExV)*100,4)
	else:
		Error2 = 'nan'

	print ">>> dHVaporization =", str(Enthalpy), "+-", str(Std_Dev), "kJ/mol"
	print ">>> Experimental =", ExV, "kJ/mol"
	print ">>> Error =", Error2,"%\n"

	string2 = '\ndHvap = ' + str(Enthalpy) + ' +- ' + str(Std_Dev) + ' kJ/mol'
	string2_1 = '\nExperimental = ' + str(ExV) + ' kJ/mol'
	string2_2 = '\nError = ' + str(Error2) + '%\n'
	fout.write(string2)
	fout.write(string2_1)
	fout.write(string2_2)

	return dHvap, Error2

def calc_dielectric(ver,prefix):
	#####################################################
	# Colecting Dieletric Constant

	dipole = calc_dipole(prefix + '.gro',prefix + '.top') # Calculating dipole for topology used
	time_ps = drift_read(ver,prefix)[2]					  # Extracting simulation time

	list_Diec = []

	for i in range(0,time_ps-20000,1000):
			f8 = open("Diec/"+ prefix + ".Diec-"+str(i)+".out","r")
			line = "text"
			while True:
					line = f8.readline()
					if not line: break
					if len(line.split())>0:
							if line.split()[0] == "Epsilon":
									break

			diec = float(line.split()[2])
			list_Diec.append(diec)


	Diec = round(float(numpy.average(list_Diec[-20:])),4)   # Set to -20
	Std_Dev = round(float(numpy.std(list_Diec[-20:])),2)    # Set to -20

	if not ExDc == 'nan':
			Error3 = round(abs((ExDc - Diec)/ExDc)*100,4)
	else:
			Error3 = 'nan'

	print ">>> Dipole from topology =", str(dipole), "D"
	print ">>> Dielectric Constant =", str(Diec), "+-", str(Std_Dev)
	print ">>> Experimental =", ExDc
	print ">>> Error =", Error3,"%\n"

	out = open('Diec/Diec-Convergence.xvg', 'w')
	t = 0
	for i in list_Diec:
			string1 = str(i).rjust(10 , ' ')
			t += 10
			out.write(str(t))
			out.write(string1)
			out.write('\n')

	out.close()

	string3_0 = '\nDipole from topology = ' + str(dipole) + ' D'
	string3 = '\nDielectric Constant = ' + str(Diec) + ' +- ' + str(Std_Dev)
	string3_1 = '\nExperimental = ' + str(ExDc)
	string3_2 = '\nError = ' + str(Error3) + '%\n'
	fout.write(string3_0)
	fout.write(string3)
	fout.write(string3_1)
	fout.write(string3_2)

	return Diec, Error3

def calc_TexpC(ver,prefix):

	simulations = ['Tlow','T0', 'Thigh']

	den_values = []
	temp_values = [Tbelow, T, Tabove]

	for i in simulations:
		f0 = open("TexpC/"+ prefix + ".Den-"+str(i)+".out","r")
		den = "text"
		while True:
			den = f0.readline()
			if not den: break
			if len(den.split())>0:
				if den.split()[0] == "Density":
					break

		density = float(den.split()[1])/1000 #convert to g/cm^3
		ln = math.log(density)
		den_values.append(ln)

	slope, intercept, r_value, p_value, std_err = stats.linregress(den_values,temp_values)

	TexpC = round(-slope/1000,4) # convert from 1/K to (10^-3)/K
	Std_Dev = round(std_err/1000,2)

	if not ExE == 'nan':
		Error4 = round(abs((ExE - TexpC)/ExE)*100,4)
	else:
		Error4 = 'nan'

	#print den_values
	#print temp_values

        print ">>> Thermal Expansion Coefficient =", str(TexpC), "+-", str(Std_Dev)
        print ">>> Experimental =", ExE
        print ">>> Error =", Error4,"%\n"

	string4 = '\nThermal Expan. Coeff = ' + str(TexpC) + ' +- ' + str(Std_Dev) + ' 10^-3/K'
	string4_1 = '\nExperimental = ' + str(ExE) + ' 10^-3/K'
	string4_2 = '\nError = ' + str(Error4) + '%\n'
	fout.write(string4)
	fout.write(string4_1)
	fout.write(string4_2)

	return TexpC, Error4

def calc_Cv(ver,prefix):

	simulations = ['Tlow','T0', 'Thigh']

	energy_values = []
	temp_values = [Tbelow, T, Tabove]

	for i in simulations:
		f0 = open("Cv/"+ prefix + ".Etot-"+str(i)+".out","r")
		string = "text"
		while True:
			string = f0.readline()
			if not string: break
			if len(string.split())>0:
				if string.split()[0] == "Total":
					Energy = float(string.split()[2])
					break

		Energy = float(Energy*1000) #convert to J/mol
		energy_values.append(Energy)

	slope, intercept, r_value, p_value, std_err = stats.linregress(temp_values,energy_values)


	Cv = round(slope + Cvib, 4)
	Std_Dev = round(std_err,4)

	if not ExCv == 'nan':
		Error5 = round(abs((ExCv - Cv)/ExCv)*100,4)
	else:
		Error5 = 'nan'

	#print temp_values
	#print energy_values

        print ">>> Heat Capacity =", str(Cv), "+-", str(Std_Dev), "J/mol*K"
        print ">>> Experimental =", ExCv, "J/mol*K"
        print ">>> Error =", Error5,"%\n"


	string5 = '\nHeat Capacity = ' + str(Cv) + ' +- ' + str(Std_Dev) + ' J/mol*K'
	string5_1 = '\nExperimental = ' + str(ExCv) + ' J/mol*K'
	string5_2 = '\nError = ' + str(Error5) + '%\n'
	fout.write(string5)
	fout.write(string5_1)
	fout.write(string5_2)

	return Cv, Error5

def calc_Comp(ver,prefix):

	simulations = ['Plow','P0', 'Phigh']

	den_values = []
	press_values = []
	ln_values =[]

	for i in simulations:
		f0 = open("Comp/"+ prefix + ".Press-"+str(i)+".out","r")
		p = "text"
		while True:
			p = f0.readline()
			if not p: break
			if len(p.split())>0:
				if p.split()[0] == "Pressure":
					break


		pressure = float(p.split()[1])/(10**4) # convert to GPa (1bar = 0.0001 GPa)
		press_values.append(pressure)

	for i in simulations:
		f0 = open(prefix + "."+str(i)+".gro","r")
		den = "text"
		den = f0.readlines()[-1]
		line = den.strip('\n')
		box_vector = float(line.split()[0])

		vol = box_vector**3

		mass = (MW*N)/((6.02)*0.1)

		density = mass/vol
		ln = math.log(density)
		den_values.append(density)
		ln_values.append(ln)


	slope, intercept, r_value, p_value, std_err = stats.linregress(press_values,ln_values)

	Comp = round(slope,4)
	Std_Dev = round(std_err,4)

	if not ExC == 'nan':
		Error6 = round(abs((ExC - Comp)/ExC)*100,4)
	else:
		Error6 = 'xabu'

	#print press_values, 'GPa'
	#print den_values, 'den'

        print ">>> Isothermal Compressibility =", str(Comp), "+-", str(Std_Dev), "1/GPa"
        print ">>> Experimental =", ExC, "1/GPa"
        print ">>> Error =", Error6,"%\n"

	string6 = '\nIsothermal Compressibility = ' + str(Comp) + ' +- ' + str(Std_Dev) + ' 1/GPa'
	string6_1 = '\nExperimental = ' + str(ExC) + ' 1/GPa'
	string6_2 = '\nError = ' + str(Error6) + '%\n'
	fout.write(string6)
	fout.write(string6_1)
	fout.write(string6_2)

	return Comp, Error6

def main(ver,prefix):

	#os.system('rm ' + prefix + '*xvg ' + prefix + '*out Density/* Comp/* dHvap/* TexpC/* Cv/* MSD/*')
	#calc_prop(ver,prefix,block)

	value = drift_output(ver,prefix,DF,name,T)

	print "\n####################################################################\n"

	print "Set of Energies:"

	DeltaD = calc_diffusion(ver,prefix)
	Den = calc_den(ver,prefix)
	dHvap = calc_dHvap(ver,prefix,block)
	Diec = calc_dielectric(ver,prefix)
	calc_TexpC(ver,prefix)
	calc_Cv(ver,prefix)
	calc_Comp(ver,prefix)

	fout.close()

main(ver,prefix)
