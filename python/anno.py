#!/usr/bin/env python

import re
import pandas as pd
import numpy as np
import time
from info import *
import os


#This file is for defining exon boundaires on the mutations from the COSMIC databse.
def main():

    dic = getExonDict()


    #Information is loaded from report.csv and new output file is generated for furthat annotation 
    with open ('output_new_with_absolute_position.csv', 'w') as ph:

        with open (base_path+'/report.csv', 'r') as f:

            print (f.readline().strip()+",closeToExon,optimality change", file=ph)
            for index, line in enumerate(f,1):
                try:
                    fields = line.split(',')
                    ens = fields[1]
                    genomePosition = int(fields[3])
                except IndexError:
                    print (line)
                    print("-------------------------------------------------")
                    print("line = "+str(index))

                exon_boundary = checkSpliceProximity(ens, genomePosition, dic, notfound)
                
                w = fields[9]
                m = fields[10]
                try:
                    diff = optChange(w,m)
                    print (line.strip()+','+str(exon_boundary)+','+str.format('{0:.3f}',diff), file=ph)
                except KeyError:
                    print("line number +"+str(index)+"\n"+line+'\n'+fields[11])

    # Print information not found into the things th
    with open ('notfoundinannotation.txt', 'w') as f:
        for line in notfound:
            print (line, file=f)


#Return integer of closest nucleotide boundary
def checkSpliceProximity(ens_name,position, dic, notfound):

    if dic.get(ens_name):
        for positions in dic[ens_name]:
            if position >= positions[0] and positions[1] >= position:
                l = position - positions[0]
                r = positions[1] - position

                return min(l,r)

    #Otherwise put in error file
    else:
        notfound.append(ens_name)


#Returns a dictionary with the start and end positions of the exons
def getExonDict():

    dic = {}
    #Absolute path here because Homo_sapiens.GRCh38.88.gff3 is 404M
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
    
    return dic


notfound = []


base_path = os.path.dirname(os.path.abspath(__file__)).rsplit("/",2)[0]

if __name__ == "__main__":
    main()
