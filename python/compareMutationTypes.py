#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import generateMutations as gm
import scipy.stats.stats as sp
import multiprocessing as mp
import time
from info import *

def main():
    cancerMutationCounts = getMutationCountsFromFile(base_path+'/output.csv',9,10)

    probMutationCounts = getMutationCountsFromFile(base_path+'/generatedMutations/generatedSynMutationsWithProb.csv',0,1)

    ranMutationCounts = getMutationCountsFromFile(base_path+'/generatedMutations/generatedMutations.csv',0,1)
    
    compareAgainstCancerProb = compareMutationCounts(cancerMutationCounts,probMutationCounts) 
    makeScatterPlot(compareAgainstCancerProb,'/home/chris/Dropbox/BIN_3005/R_graphs/figures/probSynCancer.png')
    
    compareAgainstCancerRan = compareMutationCounts(cancerMutationCounts,ranMutationCounts)
    
    makeScatterPlot(compareAgainstCancerRan, base_path+'/R_graphs/figures/ranSynCancer.png')
    
    makeScatterForRscore(createChartData(), base_path+'/R_graphs/figures/rScoreCampairaisons.png')
    #makeScatterPlot(compareAgainstCancerSyn)    
    
    
def createSynDic():
    l = {}
    for codon in syn:
        for sin in syn[codon]:
            l[codon+':'+sin] = 0
    return l



def getMutationCountsFromFile(filepathandFile,wi,mi):
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


def getMutationCountsFromMemory(mutations):
    l = createSynDic()
    for line in mutations: 
        vals = line.split(',')
        k = vals[0]+':'+vals[1]
        l[k] += 1
    return l


def compareMutationCounts(counts1,counts2):

      
    sum1 = sum(counts1.values())
    sum2 = sum(counts2.values())
    
    bigger = max(sum1, sum2)
    smaller = min(sum1, sum2)

    factor = smaller/bigger

    adjust = None
    if bigger == sum1:
        adjust = counts1
    else:
        adjust = counts2
        
    if adjust is not None:
        for mutations in adjust:
            adjust[mutations] = int(adjust[mutations]*factor)
    
    a = [[],[],[]] 
    for k in counts1:
        a[0].append(k)
        a[1].append(counts1[k])
        a[2].append(counts2[k])

    
    return a


def run2Syn(number):
    
    seq1 = gm.getSequencesFromReference()
    mutations1 = gm.getMutations()
    report1 = gm.generateSynonymousMutations(number, seq1, mutations1)
    counts1 = getMutationCountsFromMemory(report1)
    
    seq2 = gm.getSequencesFromReference()
    mutations2 = gm.getMutations()
    report2 = gm.generateSynonymousMutations(number, seq2, mutations2)
    counts2 = getMutationCountsFromMemory(report2)

    two_counts = compareMutationCounts(counts1,counts2)

    are = sp.pearsonr(two_counts[1],two_counts[2])

    print("For "+str(number)+" mutations R value is "+str(are[0])+" p-value "+str(are[1]))
    return (number, are[0])

def run2Ran(number):
    
    seq1 = gm.getSequencesFromReference()
    mutations1 = gm.getMutations()
    report1 = gm.generateRandomMutations(number, seq1, mutations1)
    counts1 = getMutationCountsFromMemory(report1)
    
    seq2 = gm.getSequencesFromReference()
    mutations2 = gm.getMutations()
    report2 = gm.generateRandomMutations(number, seq2, mutations2)
    counts2 = getMutationCountsFromMemory(report2)

    two_counts = compareMutationCounts(counts1,counts2)

    are = sp.pearsonr(two_counts[1],two_counts[2])

    print("For "+str(number)+" mutations R value is "+str(are[0])+" p-value "+str(are[1]))
    return (number, are[0])


def createChartData():
    pool = mp.Pool(processes=8)
    amounts = [1000,10000,100000,1000000]
    
    results = []
    results1 = []
    tic = time.clock() 
    for nums in amounts:
        results.append(pool.apply_async(run2Syn,args=(nums,)))
        results1.append(pool.apply_async(run2Ran,args=(nums,)))
    
    output = [p.get() for p in results]
    output1 = [p.get() for p in results1]
    value = [[[],[]],[[],[]]]
    for thing in output:
        value[0][0].append(thing[0])
        value[0][1].append(thing[1])
    
    for thing in output1:
        value[1][0].append(thing[0])
        value[1][1].append(thing[1])
    toc = time.clock()
    return value 

def makeScatterPlot(values, filepath):
    
    x = np.fromiter(values[2],dtype=float)
    y = np.fromiter(values[1],dtype=float)

    plt.scatter(x,y)
    plt.savefig(filepath)
    plt.clf()


def makeScatterPlot2(values0, values1, filepath):

    x1 = np.fromiter(values0[2],dtype=float)
    y1 = np.fromiter(values0[1],dtype=float)

    #a = plt.scatter(x,y)

    x2 = np.fromiter(values1[2],dtype=float)
    y2 = np.fromiter(values1[1],dtype=float)
    #b = plt.scatter(x,y)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.scatter(x1, y1, s=10, c='b', marker="s", label='ts/tv model')
    ax1.scatter(x2,y2, s=10, c='r', marker="o", label='simple model')
    plt.plot( [0,30000],[0,30000] )
    plt.legend(loc='upper right')
    plt.savefig(filepath)
    plt.clf()

def makeScatterForRscore(data, filepath):

    xboth = np.fromiter(data[0],dtype=float)
    ysyn = np.fromiter(data[1],dtype=float)
    yran = np.fromiter(data[2],dtype=float)
    ti = "Sampling of Background Mutations"
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.xscale('log')

    ax1.scatter(xboth, ysyn, s=10, c='b', marker="s", label='ts/tv model')
    ax1.scatter(xboth,yran, s=10, c='r', marker="o", label='simple model')
    xvals = np.array([1000,10000,100000,1000000])
    my_xticks = ['1000','10,000','100,000','1,000,000']
    yvals = np.array([0.8,0.85,0.9,0.95,1])
    plt.xticks(xvals, my_xticks)
    plt.yticks(np.arange(0.8,1.0125,0.05))
    ax1.plot(xboth, ysyn)
    ax1.plot(xboth,yran)
    plt.title(ti)
    plt.xlabel('Number of mutations generated')
    plt.ylabel('Correlation Coefficient')
    plt.ylim(0.8,1.0125)
    plt.legend(loc='lower right');
    plt.savefig(filepath)

def writeComparisonGraphToFile(info):

    reportWriter('mutation_generation_comparison_graph_data.csv','1000,10000,100000,1000000')

base_path = os.path.dirname(os.path.abspath(__file__)).rsplit("/",2)[0]

if __name__ == "__main__":
    main()
