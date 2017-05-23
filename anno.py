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
            if abs(position - positions[0]) <= 100 or abs(positions[1] -
                                                         position) <= 100:
                print(ens_name)
                print("Value of mutation position "+str(position))
                print("Value in tuple positions[0] = "+str(positions[0])+"""
                positions[1] = """+str(positions[1]))
                return True

            else:
                return False

    else:
        notfound.append(ens_name)


with open ('output.csv', 'w') as ph:

    with open ('report.csv', 'r') as f:

        print (f.readline().strip()+",closeToExon", file=ph)
        within30 = ''
        for index, line in enumerate(f,1):
            within30 = 'no'
            try:
                fields = line.split(',')
                ens = fields[0]
                genomePosition = int(fields[2])
            except IndexError:
                print (line)
                print("-------------------------------------------------")
                print("line = "+str(index))

            if checkSpliceProximity(ens, genomePosition):
                within30 = 'yes'
            print (line.strip()+','+within30, file=ph)
            # number += 1

with open ('notfoundinannotation.txt', 'w') as f:
    for line in notfound:
        print (line, file=f)
