#!/usr/bin/env python

import re
import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt

dic = {}
notfound = []

data = pd.read_csv('output.csv', sep=',', skiprows=200,
     names=["seqid","chrom","chrom_pos","cds_pos","nuc_pos_codon","oNuc","mNuc","oCodon","mCodon","type","relPos","aa","splice"]
                  )

exons = data.loc[(data['type'] == 'transversion')]


exons['relPos'].hist()
plt.show()
print(exons)
