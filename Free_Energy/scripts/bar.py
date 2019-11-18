import os
import sys


ver = sys.argv[1]
prefix = sys.argv[2]

os.chdir("/home/marcelodepolo/Heterocycles/dG_solvation/" + prefix + "/")

os.system("mkdir BAR/")
os.system("cp LAMBDA_*/dhdl*xvg BAR/")
os.chdir("/home/marcelodepolo/Heterocycles/dG_solvation/" + prefix + "/BAR/")
os.system("gmx_" + ver + " bar -f dhdl-*xvg -o -oi -oh > bar.out")
