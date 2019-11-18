import os
import sys

if sys.argv[1] == '-r':
	os.system("rm -r LAMBDA*")
elif sys.argv[1] == '-f':
	os.system("rm rm *gro *edr *trr *mdp *itp *top \#* *out *tpr *log")
elif sys.argv[1] == '-all':
	os.system("rm -r LAMBDA* rm *gro *edr *trr *mdp *itp *top \#* *out *tpr *log")
else:
	sys.exit('\n>>> Please provide -r for folders, -f for files or -all for both')

