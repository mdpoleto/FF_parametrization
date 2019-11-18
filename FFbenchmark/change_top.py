import os
import sys


file = sys.argv[1]

atomtype = ['CH3', 'F', 'S', 'CH2', 'CH1', 'CL', 'CH1']

new_atomtype = 'OA'
new_atomname = 'O'


def change(file1):

	o = open('507.top', 'w')
	
	flag = 0

	dic = {}

	f = open(file1, 'r')
	lines = len(f.readlines())

	f.seek(0)

	line = f.readline()

	while line[:9] != "[ atoms ]":
		o.write(line)
		line = f.readline()

	for i in range(3):
		o.write(line)
		line = f.readline()

	while line[:9] != "[ bonds ]":
		if len(line) > 1:
			o.write(line[:10])
			if line.split()[1] in atomtype:
				flag += 1
				a_type = new_atomtype.rjust(7,' ')
				new_line = a_type + line[17:31]
				o.write(new_line)
				atomname = new_atomname + str(flag)
				a_name = atomname.rjust(7,' ')
				insert = a_name + line[38:]
				o.write(insert)
			else:
				o.write(line[10:])
			line = f.readline()
		else:
			o.write(line)
			line = f.readline()

	while True:
		o.write(line)
		line = f.readline()
		if not line: break


change(file)
