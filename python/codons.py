#! /usr/bin/env python3


from info import *
from compareMutationTypes import getMutationCountsFromFile
import argparse, sys
import pandas as pd
import numpy as np
import scipy.stats as sp
import generateMutations as gm



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", help="Codon of interest must be specified")
    parser.add_argument("-s", help="Subtype to be analyzed must be specified")
    parser.add_argument("-f", action="store_true", help="print to file")
    parser.add_argument("-o", action="store_true", help="Oncogenes")
    parser.add_argument("-t", help=" nCodon (normal codon) or mCodon (mutatnt codon)")

    args = parser.parse_args()
    
    info = [] 

    if args.c == None:
        stats = genCase(args.s,args.f,args.t,args.o)
        #print(stats)
    
    if args.c is not None and args.s is None:
        stats = treatCodon(args.c, df_ref, args.t)
        #print(stats)

    if args.c is not None and args.s is not None:
        df = df_ref[(df_ref.cancerType == args.s)]
        stats = treatCodon(args.c, df, args.t)
        #print(stats)


def genCase(sub, ptf, target, onc):
    df = df_ref
    out = {}

    # Subset by oncogenes if onco flag true
    if onc is not None:
        df = df[df['gene'].isin(oncos)]

    for codon in AAcodons:

        if sub is not None:
            df = df[(df.cancerType == sub)]
        out[codon] = treatCodon(codon, df, target)
        
    for item in out:
        if out[item] is not None:
            print("For Codon "+item+' or = '+str(out[item][0])+' pvalue = '+str(out[item][1]))
        else:
            print (item+' == None')
            
    if ptf is not None:
        if sub is None:
            sub = "full_data_set"
        if onc is not None:
            sub = "onco_" + sub

        filename = sub+'_'+target+'.csv'
        with open(filename, 'w') as f:
            print('codon,or,pv', file =f)
            for item in out:
                if out[item] is not None:
                    print(item+','+str(out[item][0])+','+str(out[item][1]), file =f)
                else:
                    print (item+',None,None', file =f)
    return out

def treatCodon(codon,df,target):
    #Get Data frame which contains codons of interest
    dfc = df[(df[target] == codon)]
    
    #Unique genelist for genes found in filtered data frame
    genes = set(df.gene.unique())
    trans = getTranscripts(geneDic, genes)
    if len(trans[0]) < 22:
        return None

    genM = gm.generateSynonymousMutations(100000, seq, trans)
    #Generated mutation data frame
    gdf = convertGMToDF(genM)
    #Generated Mutation codon data frame
    gdfc = gdf[(gdf[target] == codon)]
    
    #Codon counts for contingency table from cancer data set
    countsCodonOfInterest = len(dfc) 
    countsTheRest = len(df) - countsCodonOfInterest
    
    #Codon counts for contingency table from generated data set
    countsCodonOfInterestG = len(gdfc) 
    countsTheRestG = len(gdf) - countsCodonOfInterestG
    ct = np.array([[countsCodonOfInterest, countsTheRest], [countsCodonOfInterestG, countsTheRestG]])
    if 0 in ct:
        print("Hi")
        return None
    #print(ct)
    oddsratio, pvalue = sp.fisher_exact( ct )

    return [oddsratio, pvalue, ct]


def convertGMToDF(l):
    ol = []
    for line in l:
        w = line.split(',')
        ol.append([w[0],w[1],w[2]])
    np_arr = np.array(ol)
    df = pd.DataFrame(np_arr, columns=['nCodon','mCodon','OC'])
    return df

#Get counts of mutations away from codons returns dictionary with codon and number takes
#getMutationCountsFromFile as input

def getValues(mutations):
    d = {}
    for codon in codons:
        
        for value in mutations:
            c = value.split(':')[0]
            if c == codon:
                if d.get(c) == None:
                    d[c] = mutations[value] 
                else:
                    d[c] += mutations[value]
    return d    

def getdf():

    df = pd.read_csv('/home/chris/Dropbox/BIN_3005/output.csv')
    print(df.columns.values)
    df = df[(df.closeToExon == 'no')]
    return df

def getTranscripts(cList, subset):

    rv = set()
    for gene in subset:
        rv.update(cList[gene])
    return [rv]    

def getOncos():
    l = []
    with open("/home/chris/Dropbox/BIN_3005/code/supplemental/oncogenes.txt", 'r') as f:
        for line in f:
            l.append(line.strip())
        
    return l    




#Global variables
seq = gm.getSequencesFromReference()
geneDic = gm.getGeneDic()
df_ref = getdf()
stop = ['TGA', 'TAA', 'TAG']
AAcodons = [aa for aa in codons if aa not in stop] 
oncos = getOncos()



if __name__ == "__main__":
    main()
