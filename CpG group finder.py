#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import csv

# Read dataframe
ri = pd.read_csv(r'C:\Users\harla\OneDrive\Documents\University\MSc Bioinformatics\RP2\Raw Data\Original (May)\CpI.csv')
rg = pd.read_csv(r'C:\Users\harla\OneDrive\Documents\University\MSc Bioinformatics\RP2\Raw Data\Original (May)\CpG.csv')
# Package position of dinucleotides and ranges of CPIs


uniChr = ri.loc[:, "chr"].unique()
uniGen = rg.loc[:, "gene.id_1"].unique()




print(uniChr)


# In[ ]:





# In[3]:


i=0
j=0
c =0
rowsDin = []
rowsCpi = []
#iterate through all dinucleotides
while i<=(len(uniChr)-1):
    j=0
    ch = uniChr[i]
    posDi = rg.loc[rg['chr'] == ch,["start"]]
    pdi =[]  
    for index, rows in posDi.iterrows():  
        m =(rows.start) 
        pdi.append(m)
    cpiRan = ri.loc[ri['chr'] == ch, ['start', 'end']].apply(tuple, axis=1)
    while j<=(len(pdi)-1):
        di = pdi[j]
        
        row = rg.loc[(rg['chr']== ch) & (rg['start'] == di)]
        
        c = 0
        # iterate through all ranges
        for rs in cpiRan:
            st, en = rs
            # is dinucleotide position within range? if so, add to list of di's found in CPI range and go to next di
            if (st <= di & di <= en):
                rowsCpi.append(row)
                break

             #if not, add one to count and check next range
            else:
                c += 1

         #if di not in any range, add to
        
        if (c == (len(cpiRan))):
            rowsDin.append(row)
        j+=1
        
    i+= 1
    print(i)

df = pd.concat(rowsDin)
counts_gene = df.groupby("gene.id_1")["gene.id_1"].transform(len)
mask = (counts_gene > 5)
df[mask].to_csv('test.csv', encoding='utf-8', index=False)
rg = pd.read_csv('test.csv')
         


# In[42]:


import pandas as pd
import numpy
import scipy
from scipy import stats
i=0
k=0
l=1
ind = 0
rowsGroup = []
de = pd.DataFrame()
de = pd.DataFrame(columns=['chr', 'start', 'end', 'meth.diff', 'N', '% CpG', 'gene', "Z_score", "weighted P-value", "overlap_gene", "overlap_promoter"])
#iterate through all genes
while i<=(len(uniChr)-1):
    print(i)
    j=0
    ch = uniChr[i]
    uniGen = rg.loc[rg['chr'] == ch, "gene.id_1"].drop_duplicates()
    ugen = uniGen.tolist()
    if '.' in ugen: ugen.remove('.')
    while j<=(len(ugen)-1):
        gen = ugen[j]
        
        pos = rg.loc[(rg['chr'] == ch) & (rg['gene.id_1'] == gen), "start"].tolist()
        pos.sort()
        pval = rg.loc[(rg['chr'] == ch) & (rg['gene.id_1'] == gen), "qvalue"].values
        pv = scipy.stats.combine_pvalues(pval, method='stouffer', weights=None)
        meth = rg.loc[rg['gene.id_1'] == gen, "meth.diff"]
        og = rg.loc[(rg['chr'] == ch) & (rg['gene.id_1'] == gen), "overlap_gene_1"].tolist()
        op = rg.loc[(rg['chr'] == ch) & (rg['gene.id_1'] == gen), "overlap_promoter_1"].tolist()
        if (pos[-1]>pos[0]):
            a = ch
            b = pos[0]
            c = pos[-1]
            d = (sum(meth)/float(len(meth)))
            e = len(pos)
            f = (e/(((c - b)/2)+1)*100)
            g = pv[0]
            h = pv[1]
            aa = 'N'
            bb = 'N'
            if 1 in og:
                aa = 'Y'
            if 1 in op:
                bb = 'Y'
            
            if (f >=5) & (-25>=d<=25):
                de.loc[ind] = [a, b, c, d, e, f, gen, g, h, aa, bb]
                rowsGroup.append(de)
                ind += 1
        

        
        if len(pos)<=3:
            j+=1

        count  = (pos[0]-1)


        # iterate through all dinucleotides associated with gene, count nucleotides between CpGs
        for p in pos:
             
            count +=1
            c = 0
            
            while p != count:
                    count += 1
                    c += 1
                    if c == 50:
                        break
            if c == 50:
                    continue
            # is dinucleotide position within range? if so, add to list of di's found in CPI range and go to next di
            if (p == count):
                k +=1


            # if not, add one to count and check next range

        # if di not in any range, add to 

        l = count - (pos[0]-1)
        j+=1
        
        
        
    i+=1
df = pd.concat(rowsGroup)
df.to_csv('CpG_50.csv', encoding='utf-8', index=True)


# In[41]:


print(df)


# In[39]:


ogs = rg.loc[(rg['chr'] == ch) & (rg['gene.id_1'] == gen), "overlap_gene_1"].tolist()

print(ogs)
if 1 in ogs:
    aa = 'Y'

print(aa)


# In[24]:


type(og)


# In[11]:


print(op)

