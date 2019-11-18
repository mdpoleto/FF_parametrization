import os
import sys

print "|-------------------------------------------------------------------------------|"
print "|---######################### MD Technologies ###############################---|"

name = str(sys.argv[1])
code = str(sys.argv[2])
temp = str(sys.argv[3])
DF = str(sys.argv[4])
Den = str(sys.argv[5])
Vap = str(sys.argv[6])
TexpC = str(sys.argv[7])
Comp = str(sys.argv[8])
DieC = str(sys.argv[9])
Debye = str(sys.argv[10])
Cp = str(sys.argv[11])
Cv = str(sys.argv[12])
MW = str(sys.argv[13])



header = "NAME                          CODE      TEMP      DegFREE    Density    Vaporization    Th.exp.Coeff    Compressibility     DieC    Dipole    Cp      Cv        MW"
string0 = name.ljust(30, " ")
string1 = code.ljust(9, " ")
string2 = temp.ljust(13, " ")
string3 = DF.ljust(9, " ")
string4 = Den.ljust(15, " ")
string5 = Vap.ljust(16, " ")
string6 = TexpC.ljust(18, " ")
string7 = Comp.ljust(14, " ")
string8 = DieC.ljust(9, " ")
string9 = Debye.ljust(8, " ")
string10 = Cp.ljust(9, " ")
string11 = Cv.ljust(8, " ")
string12 = MW.ljust(1, " ")
paragraph = "\n"

def header_check(header, paragraph):
	if not os.path.exists("smGROMOS-prop.dab"):
		file = "smGROMOS-prop.dab"
		file_h = open(file, 'w')
		file_h.write(header)
		file_h.write(paragraph)
		file_h.close()
		print "\n>>> File smGROMOS-prop.dab created!"
		
	else:
		print "\n>>> Opening smGROMOS-prop.dab..."

def insert_content(string0, string1, string2, string3, string4, string5, string6, string7, string8, string9, string10, string11, string12):
		
	file = "smGROMOS-prop.dab"
	file_h = open(file, 'a+')
	file_h.write(string0)
	file_h.write(string1)
	file_h.write(string2)
	file_h.write(string3)
	file_h.write(string4)
	file_h.write(string5)
	file_h.write(string6)
	file_h.write(string7)
	file_h.write(string8)
	file_h.write(string9)
	file_h.write(string10)
	file_h.write(string11)
	file_h.write(string12)
	file_h.write(paragraph)
	file_h.close()
	print ">>> Properties of " +str(name)+" sucessfully included into database!\n"

header_check(header,paragraph)

if not name in open("smGROMOS-prop.dab").read():
	insert_content(string0, string1, string2, string3, string4, string5, string6, string7, string8, string9, string10, string11, string12)
else:
	print "\n>>> WARNING!!! There is already a entry with the same name as yours. Please, check it!\n"
