#!usr/bin/env python3

import random 
from collections import Counter
from checkMutations2 import makeCodon, makeMutatntCodon, getSequencesFromReference 
import time

def main():
    
    seq = getSequencesFromReference()
    mutations = getMutations()
    tic = time.clock()
    generateRandomMutations(1000000, seq, mutations)
    toc = time.clock()
    print('time to complete 1000000 = '+str(toc-tic))



def getMutations():

    cntGene = Counter()
    cntTran = Counter()
    
    tran = set() 


    with open('report.csv', 'r') as f:
        f.readline()

        for line in f:
            words = line.split(',')
            
            gene = words[0]
            transcript = words[1]+'_'+words[4]

            cntGene[gene] += 1
            cntTran[transcript] += 1
            
            tran.add(transcript)
    
    return [tran, cntTran, cntGene]    


def codonProbabilities():
    pass

def checkMutationType(wc,mc):
    
    syn = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
           'CGA':'R','CGT':'R','CGC':'R','CGG':'R','AGA':'R','AGG':'R',
           'AAT':'N','AAC':'N','GAT':'D','GAC':'D','TGT':'C','TGC':'C','CAA':'Q','CAG':'Q','GAA':'E','GAG':'E',
           'GGA':'G','GGT':'G','GGC':'G','GGG':'G','CAT':'H','CAC':'H','ATT':'I','ATC':'I','ATA':'I',
           'ATG':'-','TTA':'L','TTG':'L','CTA':'L','CTT':'L','CTC':'L','CTG':'L','AAA':'K','AAG':'K','TTT':'F','TTC':'F',
           'CCA':'P','CCT':'P','CCC':'P','CCG':'P','AGT':'S','AGC':'S','TCA':'S','TCT':'S','TCC':'S','TCG':'S',
           'ACA':'T','ACT':'T','ACC':'T','ACG':'T','TGG':'W','TAT':'Y','TAC':'Y',
           'GTA':'V','GTT':'V','GTC':'V','GTG':'V', 'TAA':'ST', 'TGA':'ST', 'TAG':'ST'
          }
    
    if syn[wc] == syn[mc]:
        return True
    else:
        return False
    



def generateRandomMutations(number, sequences, transcripts):
    
    report = []
    nucs = ['A', 'T', 'C', 'G']
    noomber = 0
    while number > 0:

        options = []
        # Choose random CDS from pool 
        t = random.sample(transcripts[0], 1)[0]

        # Choose random nucleotide in chosen sequence
        i = random.randint(0,(len(sequences[t])-1))
         
        # Chosen nucleotide stored as n.
        n = sequences.get(t)[i]
        
        # Wild type codon
        wc = makeCodon(sequences.get(t), i+1)
        
        #If Codon or sequence pulled are not mod3 or contain 'N' value nucleotides
        if len(wc)<3 or len(sequences.get(t)) % 3 != 0 or 'N' in sequences.get(t):
            continue

        for item in nucs:
            if item != n:
                options.append(item)
        
        mutantNuc = random.choice(options)

        mc = makeMutatntCodon(wc, i+1, mutantNuc)[1]

        if wc == mc:
            print('n equals '+n+' mutantNuc equals '+mutantNuc+' options equals')
            print(options)
            print('modulo i+1 = '+str((i+1) % 3))
            print('wc = '+wc+' mc = '+mc)
            input()
            
        
        if checkMutationType(wc, mc):
            report.append(wc+','+mc)
            number -= 1
        noomber += 1    

    print('number of iterations '+str(noomber))
    with open('generatedMutations.csv', 'w') as f:
        for line in report:
            print(line, file = f)


if __name__ == "__main__":
    main()
