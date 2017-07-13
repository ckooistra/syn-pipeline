#!/usr/bin/env python

import re
import pandas as pd
import numpy as np
import time

dic = {}
notfound = []

anno = pd.read_csv('../Reference_Genome/Homo_sapiens.GRCh38.88.gff3', sep='\t', skiprows=200,
                    names=["seqid","source","type","start","end","score","strand","phase","attributes"])
exons = anno.loc[(anno['type'] == 'exon'),['start','end', 'attributes']]

opt = {"AAA": 0.6077, "AAC": 0.9867, "AAG": 0.7746, "AAT": 0.4521, "ACA": 0.1658, "ACC": 0.1989, "ACG": 0.2188, "ACT": 0.2762,
                        "AGA": 0.1657, "AGC": 0.2409, "AGG": 0.3293, "AGT": 0.1246, "ATA": 0.1382, "ATC": 0.5591, "ATG": 0.5525, "ATT": 0.5666,
                        "CAA": 0.7459, "CAC": 0.3039, "CAG": 0.9845, "CAT": 0.1334, "CCA": 0.2763, "CCC": 0.2862, "CCG": 0.1989, "CCT": 0.3712,
                        "CGA": 0.1658, "CGC": 0.1392, "CGG": 0.1912, "CGT": 0.1934, "CTA": 0.1105, "CTC": 0.2586, "CTG": 0.3116, "CTT": 0.3591,
                        "GAA": 0.5249, "GAC": 0.5525, "GAG": 0.8309, "GAT": 0.2425, "GCA": 0.304, "GCC": 0.5967, "GCG": 0.2354, "GCT": 0.8287,
                        "GGA": 0.3039, "GGC": 0.4144, "GGG": 0.3459, "GGT": 0.1819, "GTA": 0.1658, "GTC": 0.2188, "GTG": 0.6055, "GTT": 0.3039,
                        "TAC": 0.4619, "TAT": 0.2217, "TCA": 0.1934, "TCC": 0.2387, "TCG": 0.1724, "TCT": 0.3315, "TGC": 0.8762, "TGG": 1,
                        "TGT": 0.4036, "TTA": 0.3039, "TTC": 0.4696, "TTG": 0.3182, "TTT": 0.2062}



for index, row in exons.iterrows():
    parentTranscript = row['attributes'].split(';')[0].split(':')[1]

    if parentTranscript not in dic:
        dic[parentTranscript] = []
        dic[parentTranscript].append((int(row['start']), int(row['end'])))

    else:
        dic[parentTranscript] += (row['start'], row['end'])
    # print(parentTranscript)


def checkSpliceProximity(ens_name,position):
    if dic.get(ens_name):
        for positions in dic[ens_name]:
            if abs(position - positions[0]) <= 100 or abs(positions[1] - position) <= 100:
                return True

            else:
                return False

    else:
        notfound.append(ens_name)



def optChange(wCod,mCod):
    diff = opt[mCod] - opt[wCod] 
    print('Diff '+str(diff))
    return diff

with open ('output.csv', 'w') as ph:

    with open ('report.csv', 'r') as f:

        print (f.readline().strip()+",closeToExon,optimality change", file=ph)
        within30 = ''
        for index, line in enumerate(f,1):
            within30 = 'no'
            try:
                fields = line.split(',')
                ens = fields[1]
                genomePosition = int(fields[3])
            except IndexError:
                print (line)
                print("-------------------------------------------------")
                print("line = "+str(index))

            if checkSpliceProximity(ens, genomePosition):
                within30 = 'yes'
            
            w = fields[9]
            m = fields[10]
            print("Wild = "+w+" Mutant ="+m)
            try:
                diff = optChange(w,m)
                print (line.strip()+','+within30+','+str.format('{0:.3f}',diff), file=ph)
            except KeyError:
                print("line number +"+str(index)+"\n"+line+'\n'+fields[11])
                #input("KeyError")

with open ('notfoundinannotation.txt', 'w') as f:
    for line in notfound:
        print (line, file=f)
