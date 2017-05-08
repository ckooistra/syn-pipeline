#/bin/python

from collections import Counter
import re


with open('CosmicMutantExportSilent.tsv', 'r') as f:

    number = 0
    matchNumber = 0
    unMatched = 0
    cont = 0
    report = []
    for line in f:
        values = line.split("\t")

        report = []
        number += 1

        if number == 1000:
            break
        if values[25] == 'y' or values[23] == '' or re.search('_', values[17]):
            continue

        ensemblCode = values[1]
        position = int(values[17][2:-3])
        normalNuc = values[17][-3]
        mutantNuc = values[17][-1]
        length = int(values[2])
        genome = values[23].split(':')[0]
        genomePosition = int(values[23].split('-')[1])


        matched = False

        with open('tmp', 'r') as ph:
            ph.readline()

            for ligne in ph:
                sequences = ligne.split(" ")
                seq = sequences[1].strip()

                if (len(seq) < length):
                    continue

                code = sequences[0].split(".", 1)[0]
                # print (sequences[0])

                if (code == ensemblCode and seq[position-1] == normalNuc):

                    matchNumber += 1
                    matched = True


        if not(matched):
            unMatched += 1
            print("unMatched position = "+str(position)+"\n ensemblCode ="+ensemblCode+"\n ")

        if number % 200 == 0:
            print(number)


with open ('res.txt', 'w') as f:

    f.write("Total = "+str(number)+'\n')
    f.write("Matched = "+ str(matchNumber)+'\n')
    f.write("Un Matched = "+ str(unMatched)+'\n')
    f.write("Match percentage = "+str(matchNumber/(matchNumber+unMatched))+'\n')
