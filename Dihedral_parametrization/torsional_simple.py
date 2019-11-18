import sys
import re

prefix = sys.argv[1]
output = prefix[:-4] + "_qm.dat"

def main(file,output):

	inp = open(prefix,'r')
	values_raw = inp.readlines()

	values_list = []

	for i in values_raw:
		i = float(i.split()[0].replace("\n",''))*2625.50184
		values_list.append(i)


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


#print "\n##################################################\n>>>>>> Searching for " + mode + " values in file " + prefix
main(prefix,output)
print "\n>>> Torsional profile calculated sucessfully!\n"
