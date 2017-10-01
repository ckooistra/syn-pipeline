#!/usr/bin/env python

import re
import pandas as pd
import numpy as np
import time
from info import *

def main():


    dic = {}
    notfound = []

    anno = pd.read_csv('/home/chris/hd1/Reference_Genome/Homo_sapiens.GRCh38.88.gff3', sep='\t', skiprows=200,
                        names=["seqid","source","type","start","end","score","strand","phase","attributes"])
    exons = anno.loc[(anno['type'] == 'exon'),['start','end', 'attributes']]


    ''' Insert exon start and stop coordinates into 
    dctionary dic '''
    for index, row in exons.iterrows():
        parentTranscript = row['attributes'].split(';')[0].split(':')[1]

        if parentTranscript not in dic:
            dic[parentTranscript] = []
            dic[parentTranscript].append((int(row['start']), int(row['end'])))

        else:
            dic[parentTranscript] += (row['start'], row['end'])
        # print(parentTranscript)


    #Information is loaded from the anno.py and is used for the annotation of exon proximity
    with open ('output_30nuc.csv', 'w') as ph:

        with open ('/home/chris/hd1/COSMIC_V80/report.csv', 'r') as f:

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

                if checkSpliceProximity(ens, genomePosition, dic, notfound):
                    within30 = 'yes'
                
                w = fields[9]
                m = fields[10]
                try:
                    diff = optChange(w,m)
                    print (line.strip()+','+within30+','+str.format('{0:.3f}',diff), file=ph)
                except KeyError:
                    print("line number +"+str(index)+"\n"+line+'\n'+fields[11])
                    #input("KeyError")
    # Print information not found into the things th
    with open ('notfoundinannotation.txt', 'w') as f:
        for line in notfound:
            print (line, file=f)


#Return true if called mutation is within 100 nucleotides of called mutation
def checkSpliceProximity(ens_name,position, dic, notfound):
    if dic.get(ens_name):
        for positions in dic[ens_name]:
            if abs(position - positions[0]) <= 30 or abs(positions[1] - position) <= 30:
                return True

            else:
                return False

    else:
        notfound.append(ens_name)


if __name__ == "__main__":
    main()
