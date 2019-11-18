import os
import sys

filename = sys.argv[1]


def freq_evaluator(filename):

    f = open(filename, 'r')
    line = ''

    while True:
        values = []

        line = f.readline()
        if not line: break
        if len(line.split()) > 0:
            if line.split()[0] == 'Frequencies':
                v1 = float(line.split()[2])
                v2 = float(line.split()[3])
                v3 = float(line.split()[4])

            else:
                pass
        else:
            pass

    print '\n>>>>>> Evaluating Frequencies:'

    flag = 0

    if v1 > 0 and v2 > 0 and v3 > 0:
        print '>>> All frequencies are fine!\n'
        flag = 1
    else:
        print '>>>>>> Negative frequencies!\n'
        flag = 0
        sys.exit()

    return flag

def extract_Cvib(filename):

    flag = freq_evaluator(filename)

    if flag == 1:
        f = open(filename, 'r')
        f.seek(0)
        line = ''

        while True:
            values = []

            line = f.readline()
            if not line: break
            if len(line.split()) > 0:
                if line.split()[0] == 'Vibrational' and line.split()[1] != 'temperatures:':
                    Cvib = float(line.split()[2]) # cal/(mol*K)
                    Cvib = Cvib*(4.18) # transforming to J/(mol*K)

    print ">>> Vibrational Contribution of = " + str(Cvib) + " J/(mol*K)"
    return Cvib

extract_Cvib(filename)
