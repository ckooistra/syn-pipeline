#/bin/python

from collections import Counter
import re

def makeCodon(seq, p):
    if position % 3 == 0:
        return seq[p-2:p+1]
    if position % 3 == 1:
        return seq[p:p+3]
    if position % 3 == 2:
        return seq[p-1:p+2]

def makeMutatntCodon(codon, pos, mutation):
    muantCodon = ""
    for x in range(len(codon)):
        if x == pos-1:
            muantCodon += mutation
        else:
            muantCodon += codon[x]
    return [codon, muantCodon]





sequences = {}

# Get files into memory
with open('tmp', 'r') as ph:
    ph.readline()

    for ligne in ph:
        s = ligne.split(" ")
        ec = s[0].strip()
        seq = s[1].strip()

        sequences[ec] = seq

with open('CosmicMutantExportSilent.tsv', 'r') as f:

    number = 0
    matchNumber = 0
    unMatched = 0
    report = []
    skipped = []


    for line in f:
        if number is 1000:
            break
        values = line.split("\t")

        report = []
        number += 1

        if values[25] == 'y' or re.search('_', values[17]):

            skipped.append(values[1])
            continue

        ensemblCode = values[1]
        position = int(values[17][2:-3])
        normalNuc = values[17][-3]
        mutantNuc = values[17][-1]
        length = int(values[2])


        for code in sequences:
            if ensemblCode in code and len(sequences[code]) == length:
                codons = makeMutatntCodon(makeCodon(sequences[code], position),position % 3, mutantNuc)
                report.append(ensemblCode+","+str(position)+","+str(position%3)+","+ normalNuc +"," + mutantNuc + "," + codons[0]+ "," + codons[1])

    if number % 200 == 0:
        print(number)




with open ('report.csv', 'w') as fil:
    fil.write('ENSEMBL,Position,Nuc,mNuc,Codon,Codon Position,AA\n')
    for line in report:
        fil.write(line)

with open ('skippedReport.txt', 'w') as fil:
    for line in skipped:
        fil.write(line+"\n")
