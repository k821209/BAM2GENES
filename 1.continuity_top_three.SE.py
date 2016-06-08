#!/usr/bin/python

from __future__ import print_function
import subprocess
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
sys.path.append('./')
import pysam
import kang
import math
import re 
file_bam = sys.argv[1] #'intron3000.merge.sorted.bam'
file_fa  = sys.argv[2] #'Creinhardtii_281_v5.0.fa'
dicHD2seq = kang.Fasta2dic(file_fa)

def cigar_parse(cigar):
    match = re.findall(r'(\d+)(\w)', cigar)
    return match

def get_block(array,depth_cut=0):
    lim_len_block = 100
    #depth_cut     = 0 # ... 10. ... .. .. ..
    block_list = []
    #print(len(np.shape(array)))
    if len(np.shape(array)) == 1:
        rows = 1
        block = []
        for n,j in enumerate(array):
            if j > depth_cut:
                block.append(n)
            else:
                if len(block) > lim_len_block:
                    block_list.append([block[0],block[-1]])
                    block = []
                else:
                    block = []
        if block != []:
            block_list.append([block[0],block[-1]])
    else: 
        rows, columns = np.shape(array)       
        for i in range(rows):
            earray = array[i]
            block = []
            for n,j in enumerate(earray):
                if j > depth_cut:
                    block.append(n)
                else:
                    if len(block) > lim_len_block:
                        block_list.append([i,block[0],block[-1]])
                        block = []
                    else:
                        block = []
    return block_list


rows              = len(dicHD2seq.keys())
chromosomes       = dicHD2seq.keys()
chromosomes.sort()
dicN2chr          = dict(enumerate(chromosomes))
dicChr2N          = {b:a for a,b in dicN2chr.iteritems()}
columns           = max([len(x) for x in dicHD2seq.values()])-1
continuity_matrix = np.zeros([rows,columns],dtype=np.int)
#Outfile = open('chromosome.map.txt','w')
#for a,b in dicChr2N.iteritems():
#    print(a,b,sep='\t',file=Outfile)




print('start loop')

samfile = pysam.Samfile( file_bam, "rb" )
it      = samfile.fetch()
for line in tqdm(it):#$open('temp.sam.cut'): # should be changed to zero base map
    # Check qual

    #if line.is_proper_pair == False:
    #    continue
    if line.is_duplicate   == True:
        continue
    if line.is_qcfail      == True:
        continue
    if line.is_secondary   == True:
        continue

    # Check qual end
    chromosome   = line.reference_name
    startpos     = line.reference_start # zero based
    fragmentsize = line.tlen
    qname        = line.qname
    echr         = dicChr2N[chromosome]
    cigar        = line.cigarstring
    cigarM       = cigar_parse(cigar)
    cigarstrings = [x[1] for x in cigarM]
    cigarvalues  = [x[0] for x in cigarM]
    length = 0 
    for n, cigarstring in enumerate(cigarstrings):
        if  cigarstring == 'M':
            length += int(cigarvalues[n])
        elif cigarstring == 'I':
            length     += int(cigarvalues[n])
        elif cigarstring == 'N':
            length     += int(cigarvalues[n])
    #if fragmentsize > 200:
    #    pass
    #else: continue
    endpos         = startpos + length  # minus 1 for continuity value -> removed for this time just coverage purpose

    continuity_matrix[echr,startpos:endpos] += 1  # list characteristic can utillize fragment size itself.

array_contiguity = continuity_matrix

file_pk = '/ref/analysis/pipelines/pandas_df/Creinhardtii_281_v5.5.gene.gff3.pandas.df.pk'

df_gff_cre = pd.read_pickle(file_pk)
dic = {'mRNA'       : [],
       'length'     : [],
       'total.depth': [],
       'ratio.depth': [],
       'coverage (1x)'   : [],
       'coverage (10x)'   : [],
       'coverage (30x)'   : [],
      }
genelist = set([x for x,y in df_gff_cre.index])
for genename in tqdm(genelist):
    try:
        if math.isnan(float(genename)):
            continue
    except ValueError:
        pass
    #print type(genename)
    df      = df_gff_cre.loc[genename]
    mask    = (df[2]=='CDS')
    df_mRNA = df[mask].loc['1']
    try:
        chromosome = df_mRNA[0].values[0]
    except AttributeError:
        chromosome = df_mRNA[0]
    array      = df_mRNA[[3,4]].values
    try:
        r,c        = np.shape(array)
        if c != 2 :
            print('?!')
            exit()
        covered_array = []
        for i in range(r):
         
            left       = array[i,:][0] #int(df_mRNA[3])
            right      = array[i,:][1] #int(df_mRNA[4])
            echr       = dicChr2N[chromosome]
            contiguity = list(array_contiguity[echr][left-1:right]) # continuity value require minus 1 from right pos
            covered_array += contiguity
    except ValueError:
        left  = array[0]
        right = array[1]
        echr       = dicChr2N[chromosome]
        covered_array = list(array_contiguity[echr][left-1:right]) # continuity value require minus 1 from right pos
    covered_array = np.array(covered_array)
    #covered_array_10 = covered_array-9
    #covered_array_30 = covered_array-29
    length = len(covered_array)
    if 1:
        dic['mRNA'].append(genename)
        dic['length'].append(length)
        dic['total.depth'].append(sum(covered_array))
        dic['ratio.depth'].append(float(sum(covered_array))/float(length))
        dic['coverage (1x)'].append(len((covered_array >= 1).nonzero()[0])/float(length))
        dic['coverage (10x)'].append(len((covered_array >= 10).nonzero()[0])/float(length))
        dic['coverage (30x)'].append(len((covered_array >= 30).nonzero()[0])/float(length))

df_cont = pd.DataFrame(dic)
#df_cont_ix = df_cont.set_index('mRNA')



array = df_cont.sort_values(by='total.depth',ascending=False).head(3)[['mRNA','coverage (1x)','coverage (10x)','coverage (30x)','ratio.depth','total.depth']].values.ravel()
Outfile = open(file_bam+'.topthree.txt','w')
print(file_bam, '\t'.join(map(str,array)),sep='\t',file=Outfile)
