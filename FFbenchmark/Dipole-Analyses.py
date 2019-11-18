########################################################################
# Loading libraries
import numpy
import sys

arq = sys.argv[1]

gro = str(arq) + '.gro'
top = str(arq) + '.top'

# Setting parameters
eA = 4.80320440079 # Debye

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

	module = ((dipoleX**2) + (dipoleY**2) + (dipoleZ**2))**0.5

	print "###################################################\n"

	print '>>>>>> Componente vetorial em X = '+ str(round(dipoleX,2)) + ' D'
	print '>>>>>> Componente vetorial em Y = '+ str(round(dipoleY,2)) + ' D'
	print '>>>>>> Componente vetorial em Z = '+ str(round(dipoleZ,2)) + ' D'

	print '\n>>>>>> Modulo resultante = ' + str(round(module,2)) + ' D'

dic1 = calc_dipole(gro,top)

