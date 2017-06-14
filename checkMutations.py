#/bin/python

from collections import Counter
import re

def makeCodon(seq, p):
    if p% 3 == 0:
        return seq[p-2:p+1]
    if p% 3 == 1:
        return seq[p:p+3]
    if p% 3 == 2:
        return seq[p-1:p+2]

def makeMutatntCodon(codon, pos, mutation):
    muantCodon = ""
    for x in range(len(codon)):
        if x == 3-(pos+1):
            muantCodon += mutation
        else:
            muantCodon += codon[x]
    return [codon, muantCodon]

def transitionOrTransversion(normalNuc, mutantNuc):
    purines = ["A", "G"]
    pyrimidines = ["T", "C"]

    if (normalNuc in purines and mutantNuc in purines) or (normalNuc in purines and mutantNuc in pyrimidines):
        return "transition"
    else:
        return "transversion"

def positionOfMutation(pos, seqLength):
   return pos/seqLength

# Dictionary to keep sequence data
sequences = {}

# Get referemce CDS files into memory as dictionary
with open('tmp', 'r') as ph:
    ph.readline()

    for ligne in ph:
        s = ligne.split(" ")
        seq = s[1].strip()
        ec = s[0].strip().split(".")[0]+"_"+str(len(seq))

        sequences[ec] = seq

with open('CosmicMutantExportSilent.tsv', 'r') as f:

    number = 0
    matchNumber = 0
    unMatched = 0
    report = []
    skipped = []
    cnt = Counter()
    
    n = 0
    for line in f:
        # if number == 1000:
        #     break

        values = line.split("\t")

        number += 1

        if values[25] == 'y' or re.search('_', values[17]) or values[17] == 'c.?' or re.search('\*', values[18]):

            skipped.append(values[1]+", not valid")
            unMatched += 1
            continue

        try:
            gene = values[0].strip()
            ensemblCode = values[1]
            position = int(values[17][2:-3])-1
            normalNuc = values[17][-3]
            mutantNuc = values[17][-1]
            length = int(values[2])
            genome = int(values[23].split(":")[0])
            genomePosition = int(values[23].split("-")[1])
            aa = values[18][-1]

        except ValueError:
            skipped.append(ensemblCode+",missing genome position")
            unMatched += 1
            continue

        code = ensemblCode+"_"+str(length)
        if sequences.get(code):
            if len(sequences.get(code))%3 != 0:
                print("Ensembl code = "+code)
                print("First 5 nucs = "+sequences[code][0:10])
                print("Last 5 nucs = "+sequences[code][-10:])
                n+=1
            codons = makeMutatntCodon(makeCodon(sequences[code],
                                                position),position+1 % 3, mutantNuc)
            report.append(ensemblCode+","+str(genome)+","+str(genomePosition)+","+str(position)+","+
            str(3-((position+1)%3))+","+ normalNuc +"," + mutantNuc + "," +
            codons[0]+ "," +
            codons[1]+","+transitionOrTransversion(normalNuc,mutantNuc)+","+str(round(positionOfMutation(position,
                length),3))+","+aa+"\n")

            matchNumber += 1
            cnt[gene] += 1
            
            if aa == '*' and len(sequences.get(code))%3 is not 0:
                print("Code is "+code+" modulo 3 "+str(len(sequences.get(code))%3))
                
        else:
            s = ensemblCode +' '+ str(length)+' ' 
            keys = [key for key, value in sequences.items() if ensemblCode in
             key.upper()]
            for k in keys:
                s += '---'+k+' '+str(len(sequences[k]))

            skipped.append(s)

        if number % 200 == 0:
            print(number)


print("n(number of items that did not work out to modulo 0 but the coding\
      sequence was not there)= "+str(n))
print ("Matched number = "+str(matchNumber))
print ("unMatched number = "+str(unMatched))
print ("Match percentage = "+str((matchNumber/(unMatched+matchNumber)*100)))


with open ('report.csv', 'w') as fil:
    fil.write('ENSEMBL,Genome,Genome Position,CDS Position,nucOfCodon,nuc,mNuc,nCodon,mCodontion,TransitionOrTransversion,position in AA chain %,AA\n')
    for line in report:
        fil.write(line)

with open ('skippedReport.txt', 'w') as fil:
    for line in skipped:
        fil.write(line+"\n")

with open ("geneCount.txt", 'w') as fil:
    for key in cnt.most_common():
        fil.write(key[0]+" "+str(key[1])+'\n')
