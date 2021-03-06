#/bin/python

from collections import Counter
import re,os

def main():
    
    s = getSequencesFromReference()
    getMutationInfo(s)

#Return codon from position given as parameter of CDS
def makeCodon(seq, p):
    if p% 3 == 0:
        return seq[p-3:p]
    if p% 3 == 1:
        return seq[p-1:p+2]
    if p% 3 == 2:
        return seq[p-2:p+1]

#Make mutant codon based on mutation position 
def makeMutatntCodon(codon, pos, mutation):

    pos = pos % 3
    muantCodon = ""
    
    if pos is 0:
        return codon[0:2]+mutation
    if pos is 1:
        return mutation+codon[1:]
    if pos is 2:
        return codon[0]+mutation+codon[1]

#Given two codons returns transition or transversion
def transitionOrTransversion(normalNuc, mutantNuc):
    purines = ["A", "G"]
    pyrimidines = ["T", "C"]

    if (normalNuc in purines and mutantNuc in purines) or (normalNuc in pyrimidines and mutantNuc in pyrimidines):
        return "transition"
    else:
        return "transversion"

#Reutrns relative position of mutation position on the CDS
def positionOfMutation(pos, seqLength):
   return pos/seqLength

# Get sequences and return dictionary with python CDS sequences
# Sequence ENSEMBL ID and length....ENSxxxxx_length
def getSequencesFromReference():
    sequences = {}
    
    with open(base_path+'/referenceGenomeCDS.txt', 'r') as f:
        f.readline()

        for line in f:
            
            s = line.split(" ")
            seq = s[1].strip()
            ec = s[0].strip().split(".")[0]+"_"+str(len(seq))

            sequences[ec] = seq
    
    return sequences

def reportWriter(name, title, info):
    with open(name, 'w') as f:
        print(title, file=f)
        for line in info:
            print(line, file=f)


#Filters mutations by checking that they can be mapped.  If for some it can't it will print out the mutation
#to one of the error files,  this return report.csv which then contains the filtered results
def getMutationInfo(sequences):
    
    #This is where you put your copy of the COSMIC v80 database that has been filtered for mutations.    
    with open('CosmicMutantExportSilent.tsv', 'r') as f:

        number = 0
        matchNumber = 0
        unMatched = 0
        report = []
        skipped = []
        cnt = Counter()
        mod = []

        n = 0
        for line in f:

            values = line.split("\t")

            number += 1
            # If known snip, or if insertion, stop codon or if sequence change not known skip
            if values[25] == 'y' or re.search('_', values[17]) or values[17] == 'c.?' or re.search('\*|p.*X', values[18]):

                skipped.append(values[1]+", not valid")
                unMatched += 1
                continue

            try:
                gene = values[0].strip()
                ensemblCode = values[1]
                position = int(values[17][2:-3]) #Position of the mutation
                normalNuc = values[17][-3]
                mutantNuc = values[17][-1]
                length = int(values[2])
                genome = int(values[23].split(":")[0])
                genomePosition = int(values[23].split("-")[1])
                aa = values[18][-1]
                tip = values[7]
                primary_histol = values[11]

            except ValueError:
                skipped.append(ensemblCode+",missing genome position")
                unMatched += 1
                continue

            code = ensemblCode+"_"+str(length)
            if sequences.get(code):
                if len(sequences.get(code))%3 != 0:
                    mod.append("Ensembl code = "+code+"First 10 nucs = "+sequences[code][0:10]+"+Last 5 nucs = "+sequences[code][-10:])

                codons = []
                codons.append(makeCodon(sequences[code], position))
                codons.append(makeMutatntCodon(codons[0], position, mutantNuc))
                report.append(gene+','+ensemblCode+","+str(genome)+","+str(genomePosition)+
                ","+str(len(sequences[code]))+","+str(position)+","+str(3-((position)%3))+
                ","+ normalNuc +"," + mutantNuc + "," + codons[0]+ ","+codons[1]+","+
                transitionOrTransversion(normalNuc,mutantNuc)+","+
                str(round(positionOfMutation(position, length),3))+","+aa+','+tip+','+primary_histol+"\n")


                matchNumber += 1
                cnt[gene] += 1

                if aa == '*' and len(sequences.get(code))%3 is not 0:
                    print("Code is "+code+" modulo 3 "+str(len(sequences.get(code))%3))

            #If the code is not found then we put it into the skipped report
            else:
                pass
                s = ensemblCode +' '+ str(length)+' '
                #keys = [key for key, value in sequences.items() if ensemblCode in
                #key.upper()]
                #for k in keys:
                 #s += '---'+k+' '+str(len(sequences[k]))

                skipped.append(s)

            if number % 200 == 0:
                print(number)


    print("n(number of items that did not work out to modulo 0 but the coding\
          sequence was not there)= "+str(n))
    print ("Matched number = "+str(matchNumber))
    print ("unMatched number = "+str(unMatched))
    print ("Match percentage = "+str((matchNumber/(unMatched+matchNumber)*100)))


    with open ('report.csv', 'w') as fil:
        fil.write('gene,ENSEMBL,Genome,Genome Position,LengthTranscript,CDS Position,nucOfCodon,nuc,mNuc,nCodon,mCodontion,TransitionOrTransversion,position_AA,AA,cancerType,histol\n')
        for line in report:
            fil.write(line)

    with open ('skippedReport.txt', 'w') as fil:
        for line in skipped:
            fil.write(line+"\n")

    with open ("geneCount.txt", 'w') as fil:
        for key in cnt.most_common():
            fil.write(key[0]+" "+str(key[1])+'\n')

    with open ("notMod.txt", 'w') as fil:
            for line in mod:
                        fil.write(line+'\n')

base_path = os.path.dirname(os.path.abspath(__file__)).rsplit("/",2)[0]

if __name__ == "__main__":
    main()
