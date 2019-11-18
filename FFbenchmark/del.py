import os
import sys

if sys.argv[1] == '-r':
	os.system("rm -r dHvap/ Density/ Comp/ Cv/ MSD/ TexpC/ Diec/")
elif sys.argv[1] == '-f':
	os.system("rm *xvg *gro *top *itp *mdp *log *tpr *edr *trr *xtc *log *cpt *dat *done *out \#*")
elif sys.argv[1] == '-all':
	os.system("rm -r dHvap/ Density/ Comp/ Cv/ MSD/ TexpC/ Diec/ *xvg *gro *top *itp *mdp *log *tpr *edr *trr *xtc *log *cpt *dat *done *out \#*")
else:
	sys.exit('\n>>> Please provide -r for folders, -f for files or -all for both')

