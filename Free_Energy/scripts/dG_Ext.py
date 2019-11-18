import os
import sys

ver = sys.argv[1]
prefix = sys.argv[2]
lamb = sys.argv[3]


if os.path.exists(prefix + ".SIM.gro"):
    	convertcommand = "gmx_" + ver + " convert-tpr -s " + prefix + ".SIM.tpr -until 3000 -o " + prefix + ".SIM.tpr"
	os.system(convertcommand)

	min1command = "gmx_" + ver + " mdrun -deffnm " + prefix + ".SIM -dhdl dhdl-"+ str(lamb)+".xvg -nb gpu -cpi " + prefix + ".SIM.cpt -append"
        os.system(min1command)

print "\n>>> SIM completed!\n"



