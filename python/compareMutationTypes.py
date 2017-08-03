#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from info import *


def main():
    pass    


def createSynDic():
    l = {}
    for codon in syn:
        for sin in syn[codon]:
            l[codon+':'+sin] = 0
    return l


def getMutationCounts(mutations):
    l = createSynDic()

def getMutations(filepathandFile,wi,mi):
    l = createSynDic()
    count = 0
    with open(filepathandFile, 'r') as f:
        f.readline()
        for num, line in enumerate(f):
            vals = line.split(',')
            if vals[wi] == vals[mi] or synAA[vals[wi]] != synAA[vals[mi]]:
                #input ('Line # '+str(num)+' '+vals[9]+' '+vals[10])
                count += 1
                continue
            k = vals[wi]+':'+vals[mi]
            l[k] += 1
        print(str(count))
        return l    

def compareMutations(cfilepathandFile,cwi,cmi,rfilepathandFile,rwi,rmi):
    
    cancerMutationCounts = getMutations(cfilepathandFile,cwi,cmi)
    randomMutationCounts = getMutations(rfilepathandFile,rwi,rmi)
    
    cancer_sum = sum(cancerMutationCounts.values())
    generated_sum = sum(randomMutationCounts.values())
    
    bigger = max(cancer_sum, generated_sum)
    smaller = min(cancer_sum, generated_sum)

    factor = smaller/bigger

    adjust = None
    if bigger == cancer_sum:
        adjust = cancerMutationCounts
    else:
        adjust = randomMutationCounts
        
    for mutations in adjust:
        adjust[mutations] = int(adjust[mutations]*factor)
    
    a = [[],[],[]] 
    for k in cancerMutationCounts:
        a[0].append(k)
        a[1].append(cancerMutationCounts[k])
        a[2].append(randomMutationCounts[k])

    
    return a


def makeScatterPlot(values, filepath):
    
    x = np.fromiter(values[2],dtype=float)
    y = np.fromiter(values[1],dtype=float)

    plt.scatter(x,y)
    plt.show()



if __name__ == "__main__":
    main()
        
