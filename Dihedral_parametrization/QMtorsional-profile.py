import sys
import re

#############################################################
# This script was writen to create a graph for torsional
# profile based on Gaussian outputs. It works with HF and MP2
# calculations type (defined on line 17).
#
#   Usage: python QMtorsional-profile.py MODE gaussian_output.out
# Example: python QMtorsional-profile.py MP2 ben.out
#
#############################################################

mode = sys.argv[1]
prefix = sys.argv[2]

######################################################
output = prefix[:-4]+"_qm.dat"   # defining output name

######################################################

if mode == 'HF' or mode == 'MP2':
	pass
else:
	sys.exit("\n>>> Please select the first argument as HF or MP2.")


def check(file):
	f = open(file,"r")
	string = "text"

	values_raw = []

	flag = True
	flag_2 = False

	while flag == True:
		string = f.readline()
		if not string:break
		if string.split() > 0:

			if mode == 'MP2':
				if "\MP2=" in string:
					flag_2 = True
					values_raw.append(string)

				elif "\RMSD=" in string:
					values_raw.append(string)
					flag_2 = False
					break

				else:
					pass

				if flag_2 == True:
					if not string in values_raw:
						values_raw.append(string)
				else:
					pass

			elif mode == 'HF':
				if "\HF=" in string:
					flag_2 = True
					values_raw.append(string)

				elif "RMSD=" in string:
					values_raw.append(string)
					flag_2 = False
					break

				else:
					pass

				if flag_2 == True:
					if not string in values_raw:
						values_raw.append(string)
				else:
					pass
			else:
				pass


	if mode == 'MP2':
		var = ''.join(values_raw)
		var = var.replace('\n', ' ').replace('\r', '').replace(' ','')
		v1 = re.search("\MP2=(.+)\RMSD",var)
	elif mode == 'HF':
		var = ''.join(values_raw)
		var = var.replace('\n', ' ').replace('\r', '').replace(' ','')
		v1 = re.search("\HF=(.+)RMSD",var)
	else:
		pass

	print v1.groups(1)[0]

	values_list_HF = v1.groups(1)[0].split(',') # in Hartree
	values_list = []

	for i in values_list_HF:                    # converting to kJ/mol
		v = float(i)*2625.50184
		values_list.append(v)

	return values_list

def main(file,output):

	values_list = check(file)
	minor = min(values_list)

	values = {}
	angle = 0

	step = 360/(len(values_list)-1)

	for i in values_list:
		energy = i - minor
		values[angle] = "{:10.4f}".format(energy)
		angle += step


	f = open(output,"w")

	for a in sorted(values.keys()):
		angle = str(a).ljust(3," ")
		energy = str(values[a])

		f.write(angle)
		f.write(energy)
		f.write('\n')


print "\n##################################################\n>>>>>> Searching for " + mode + " values in file " + prefix
main(prefix,output)
print "\n>>> Torsional profile calculated sucessfully!\n"
