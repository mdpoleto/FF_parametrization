# USAGE: 'pymol -qc rotate.py prefix.pdb A1 A2 A3 A4'
# beeing A[1-4] the atoms of the dihedral

# NOTE: some errors may appear. If the files are generated, then you are giant.

import pymol
import sys


prefix = sys.argv[3]
dt = 30 # degrees step


D1 = "name " + sys.argv[4]
D2 = "name " + sys.argv[5]
D3 = "name " + sys.argv[6]
D4 = "name " + sys.argv[7]

cmd.load(prefix)
prefix[:3]

for i in range(0,360+dt,dt):
    dihedral = "name " + D1 + ", name " + D2 + ", name " + D3 + ", name " + D4 + ", " + str(i) +""
    name = prefix[:3] + "_" + str(i) + ".pdb"

    cmd.set_dihedral(D1, D2, D3, D4, i)

    cmd.save(name)
    cmd.quit
