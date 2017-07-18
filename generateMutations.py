#!usr/bin/env python

import random 
import numpy as np
from collections import Counter
from checkMutations2 import makeCodon, makeMutatntCodon, getSequencesFromReference, transitionOrTransversion, reportWriter 
from anno import opt, optChange
import time

def main():
    
    seq = getSequencesFromReference()
    mutations = getMutations()
    tic = time.clock()
    generateRandomMutations(1000000, seq, mutations)
    toc = time.clock()
    print('time to complete 1000000 = '+str((toc-tic)/60))


synAA = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
       'CGA':'R','CGT':'R','CGC':'R','CGG':'R','AGA':'R','AGG':'R',
       'AAT':'N','AAC':'N','GAT':'D','GAC':'D','TGT':'C','TGC':'C','CAA':'Q','CAG':'Q','GAA':'E','GAG':'E',
       'GGA':'G','GGT':'G','GGC':'G','GGG':'G','CAT':'H','CAC':'H','ATT':'I','ATC':'I','ATA':'I',
       'ATG':'-','TTA':'L','TTG':'L','CTA':'L','CTT':'L','CTC':'L','CTG':'L','AAA':'K','AAG':'K','TTT':'F','TTC':'F',
       'CCA':'P','CCT':'P','CCC':'P','CCG':'P','AGT':'S','AGC':'S','TCA':'S','TCT':'S','TCC':'S','TCG':'S',
       'ACA':'T','ACT':'T','ACC':'T','ACG':'T','TGG':'W','TAT':'Y','TAC':'Y',
       'GTA':'V','GTT':'V','GTC':'V','GTG':'V', 'TAA':'ST', 'TGA':'ST', 'TAG':'ST'
      }

synGroup = { 'A':['GCT', 'GCC', 'GCA', 'GCG'], 'R':['CGA','CGT','CGC','CGG','AGA','AGG'],
            'N':['AAT','AAC'], 'D':['GAT','GAC'], 'C':['TGT','TGC'], 'Q':['CAA','CAG'], 'E':['GAA','GAG'],
            'G':['GGA','GGT','GGC','GGG'], 'H':['CAT','CAC'], 'I':['ATT','ATC','ATA'], 'M':['ATG'],
            'L':['TTA','TTG','CTA','CTT','CTC','CTG'], 'K':['AAA','AAG'], 'F':['TTT','TTC'],
            'P':['CCA','CCT','CCC','CCG'], 'S':['AGT','AGC','TCA','TCT','TCC','TCG'],
            'T':['ACA','ACT','ACC','ACG'], 'W':['TGG'], 'Y':['TAT','TAC'], 'V':['GTA','GTT','GTC','GTG'],
            'ST':['TAA','TAG','TGA'] 
           }


codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']


def createCodonLists():
    l = []
    for i in range(64):
        l.append([0] * 64)
    
    return l
#Returns a dictionary with the keys as integer positions
#and the values as a list with [wildTypeNuc, mutantNuc, transitionOrTransversion]
def checkCodonDifferences(c1, c2):
    
    rep = {}
    l1 , l2 = list(c1), list(c2)

    for i in range(3):
        if l1[i] != l2[i]:
            rep[i] = [l1[i], l2[i], transitionOrTransversion(l1[i], l2[i])]
    return rep

# Normalizes the nonzero values of a list to 1
def listNormalizer(l):

    a = []
    
    for row in l:
        newArr = []

        for elem in row:
            if sum(row)!= 0:
                newArr.append(elem/sum(row))
            else:
                newArr.append(elem)
        a.append(newArr)
    
    return a 

#Creates a numpy 2D array to provide probabilities for likely synonymous mutations
def createMutationProbMatrix(transi, transv):
    
    lists = createCodonLists()
    t = {'transition':transi, 'transversion':transv}
    for i in range(64):
        c1 = codons[i]
        
        for j in range(64):
            c2 = codons[j]
            differences = checkCodonDifferences(c1,c2)
            if c2 == c1 or synAA[c2] != synAA[c1]:
                lists[i][j] = 0
            
            else:
                prob = 1.0
                for key in differences:
                    prob *= t[differences[key][2]] 
                lists[i][j] = prob

    lists = listNormalizer(lists)

    a  = np.array(lists)

    return (a)


#Gets cleaned mutations pulled from Cosmic database
#Returns list with set of transcripts, transcriptFrequecy count and gene count 
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

#Returns a list with two vectors, one of potential codons
#the other with probabilities for the first vector
def codonWeightedChoices(probVector, probNonZero):
    
    l =[[],[]]
    
     
    #Loop over indices in probability vector
    for codonIndex in probNonZero:
        codon = codons[codonIndex]
        probability = probVector[codonIndex]
        l[0].append(codon)
        l[1].append(probability)

    
    return [tuple(l[0]), tuple(l[1])]

 
#Returns boolean valeur of synonymity
def checkMutationType(wc,mc):
    
    
    if synAA[wc] == synAA[mc]:
        return True
    else:
        return False
    


#Method to generate random mutations, takes as parameters
#Number of mutations to perform and the dictionary of sequences from reference genome
#and set of transcritps, Returns file with generated mutations including optimality change
#This method generates random mutations but discards non synonymous ones
def generateRandomMutations(number, sequences, transcripts):
    
    report = []
    nucs = np.array(['A', 'T', 'C', 'G'])
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

        
        options = nucs != n
        
        mutantNuc = random.choice(nucs[options])
        mc = makeMutatntCodon(wc, i+1, mutantNuc)
        

        if wc == 'TGA' or wc == 'TAA' or wc == 'TAG' or mc == 'TGA' or mc == 'TAA' or mc == 'TAG':
            continue
        
        oc = optChange(wc, mc)
        tOrt = transitionOrTransversion(n, mutantNuc)
        if checkMutationType(wc, mc):
            report.append(wc+','+mc+','+tOrt+','+str.format('{0:.3f}',oc))
            number -= 1
        noomber += 1    

    print('number of iterations '+str(noomber))
    
    reportWriter('generateRandomMutations.csv','WC,MC,TransitionOrTransversion,OC',report)


''' Same as above function however all mutations are synonymous calculated based on probability
vector '''

def generateSynonymousMutations(number, sequences, transcripts):
    
    report = []
    
    prob = createMutationProbMatrix(0.3,0.1)

    while number > 0:
        # Choose random CDS from pool 
        t = random.sample(transcripts[0], 1)[0]

        # Choose random nucleotide in chosen sequence
        i = random.randint(0,(len(sequences[t])-1))
         
        # Chosen nucleotide stored as n.
        n = sequences.get(t)[i]
        
        # Wild type codon
        wc = makeCodon(sequences.get(t), i+1)

        #If Codon or sequence pulled are not mod3 or contain 'N' value nucleotides
        if len(wc)<3 or len(sequences.get(t)) % 3 != 0 or 'N' in sequences.get(t) or wc == 'ATG' or synAA[wc] == 'ST' or wc == 'TGG':
            continue

        #Get Probability Vector from Probability Matrix for wildtype codon
        probVector = prob[codons.index(wc)]
        
        #Get indices of values that are not zero
        nonzeros = np.nonzero(probVector)[0]

        #Probabilities for given WildType codon (list l with codons l[0] and probabilities l[1])
        probsForWildTypeCodon = codonWeightedChoices(probVector, nonzeros)

        mc = np.random.choice(probsForWildTypeCodon[0], p=probsForWildTypeCodon[1])
        
        oc = optChange(wc, mc)

        report.append(wc+','+mc+','+str.format('{0:.3f}',oc)) 
        number -= 1

    reportWriter('generatedSynMutationsWithProb.csv', 'WC,MC,OC', report) 

if __name__ == "__main__":
    main()
