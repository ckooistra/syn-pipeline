#! /usr/bin/env python3


syn = {
    "GCT":["GCC","GCA","GCG"],"GCC":["GCT","GCA","GCG"],"GCA":["GCT","GCC","GCG"],
    "GCG":["GCT","GCC","GCA"],"CGA":["CGT","CGC","CGG","AGA","AGG"],
    "CGT":["CGA","CGC","CGG","AGA","AGG"],"CGC":["CGA","CGT","CGG","AGA","AGG"],
    "CGG":["CGA","CGT","CGC","AGA","AGG"],"AGA":["CGA","CGT","CGC","CGG","AGG"],
    "AGG":["CGA","CGT","CGC","CGG","AGA"],"AAT":["AAC"],"AAC":["AAT"],"GAT":["GAC"],
    "GAC":["GAT"],"TGT":["TGC"],"TGC":["TGT"],"CAA":["CAG"],"CAG":["CAA"],"GAA":["GAG"],
    "GAG":["GAA"],"GGA":["GGT","GGC","GGG"],"GGT":["GGA","GGC","GGG"],"GGC":["GGA","GGT","GGG"],
    "GGG":["GGA","GGT","GGC"],"CAT":["CAC"],"CAC":["CAT"],"ATT":["ATC","ATA"],"ATC":["ATT","ATA"],
    "ATA":["ATT","ATC"],"ATG":[],"TTA":["TTG","CTA","CTT","CTC","CTG"],
    "TTG":["TTA","CTA","CTT","CTC","CTG"],"CTA":["TTA","TTG","CTT","CTC","CTG"],
    "CTT":["TTA","TTG","CTA","CTC","CTG"],"CTC":["TTA","TTG","CTA","CTT","CTG"],
    "CTG":["TTA","TTG","CTA","CTT","CTC"],"AAA":["AAG"],"AAG":["AAA"],"TTT":["TTC"],"TTC":["TTT"],
    "CCA":["CCT","CCC","CCG"],"CCT":["CCA","CCC","CCG"],"CCC":["CCA","CCT","CCG"],
    "CCG":["CCA","CCT","CCC"],"AGT":["AGC","TCA","TCT","TCC","TCG"],
    "AGC":["AGT","TCA","TCT","TCC","TCG"],"TCA":["AGT","AGC","TCT","TCC","TCG"],
    "TCT":["AGT","AGC","TCA","TCC","TCG"],"TCC":["AGT","AGC","TCA","TCT","TCG"],
    "TCG":["AGT","AGC","TCA","TCT","TCC"],"ACA":["ACT","ACC","ACG"],"ACT":["ACA","ACC","ACG"],
    "ACC":["ACA","ACT","ACG"],"ACG":["ACA","ACT","ACC"],"TGG":[],"TAT":["TAC"],"TAC":["TAT"],
    "GTA":["GTT","GTC","GTG"],"GTT":["GTA","GTC","GTG"],"GTC":["GTA","GTT","GTG"],
    "GTG":["GTA","GTT","GTC"],"TAA":["TAG","TGA"],"TAG":["TAA","TGA"],"TGA":["TAA","TAG"]

}

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

tRNA = {}
def reportWriter(name, title, info):
    with open(name, 'w') as f:
        print(title, file=f)
        for line in info:
            print(line, file=f)
