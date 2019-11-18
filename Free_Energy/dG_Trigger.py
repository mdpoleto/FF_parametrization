###############################################################################
# Usage: 'python dG_Trigger.py 507 ben'
#
# This script will automatically run hydration energy simuations. For that, you
# just need your .pdb and a force field with the necessary parameters. The scripts
# will create the simulation box, and move .gro and .top files for each LAMBDA
# folder and, within each of them, run the simulations.
#
# As general rule, just copy a .pdb and a force field into this folder and change
# the few things for compatibility.
###############################################################################

import os
import sys

ver = sys.argv[1]
prefix = sys.argv[2]

lamb = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14',
 '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30',
  '31', '32', '33', '34', '35', '36', '37', '38', '39', '40']

print
os.system("pwd")
os.system("ls")
print
os.system("python scripts/dG_Prep.py " + ver + " " + prefix)

for i in lamb:
    #####################################
    # Creating lambda directories and accessing each
    os.system("mkdir LAMBDA_"+str(i))
    os.system("cp " + prefix + ".min0.gro " + prefix + ".top LAMBDA_"+str(i)+"/")
    os.system("cp -r argromos.ff/ LAMBDA_"+str(i)+"/")
    os.chdir("LAMBDA_"+str(i))

    #####################################
    # Running simulations
    os.system("pwd")
    os.system("python ../scripts/dG_Min.py " + ver + " " + prefix + " " + str(i))
    os.system("python ../scripts/dG_Equ.py " + ver + " " + prefix + " " + str(i))
    os.system("python ../scripts/dG_Prod.py " + ver + " " + prefix + " " + str(i))

    #####################################
    # Coming back for another iteration
    os.chdir("../")
