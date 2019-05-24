# This python script is for running directout
# The parameter values should be changed everytime a different case is being processed

import sys
import os

print(os.getcwd())

def main():
    # arguments
    LEVEL = '12'
    ak = '0.05'
    BO = '200.'
    RE = '100000.'
    m = '0.1'
    B = '0'
    UstarRATIO = '0.8'
    NUMBER = 32
    TIME = 2
    # NUMBER = int(argv[1])
    # TIME = int(argv[2])
    if len(sys.argv) != 3:
        print('Parameter number wrong!\n')
    for i in range(1, TIME*NUMBER):
        Snapshot = str(i*1/NUMBER)
        cmd = './directout'+'\t'+LEVEL+'\t'+ak+'\t'+BO+'\t'+RE+'\t'+m+'\t'+B+'\t'+UstarRATIO+'\t'+Snapshot
        print(cmd)
        os.system(cmd)

if __name__ == "__main__":
    main()

print("DONE")
