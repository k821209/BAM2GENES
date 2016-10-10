#!/usr/bin/python

from __future__ import print_function
import subprocess
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
sys.path.append('/mnt/c/ubuntu.download/pipelines/')
import pysam
import kang
import math
file_bam = sys.argv[1] #'intron3000.merge.sorted.bam'
file_fa = sys.argv[2]  #'Creinhardtii_281_v5.0.fa'
file_pk = sys.argv[3] 
dicHD2seq = kang.Fasta2dic(file_fa)

def get_block(array,depth_cut=0):
    lim_len_block = 100 # size of read fragment
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
match_matrix      = np.zeros([rows,columns],dtype=np.int)
#Outfile = open('chromosome.map.txt','w')
#for a,b in dicChr2N.iteritems():
#    print(a,b,sep='\t',file=Outfile)


array_contiguity = continuity_matrix


df_gff_cre = pd.read_pickle(file_pk)
dic = {'mRNA'       : [],
       'length'     : [],
       'total.depth': [],
       'ratio.depth': [],
       'coverage (1x)'   : [],
       'coverage (10x)'   : [],
       'coverage (30x)'   : [],
       'match' : [],
       'match.ratio' :[]
      }
genelist = set([x for x,y in df_gff_cre.index])
for genename in tqdm(genelist):
    try:
        if math.isnan(float(genename)):
            continue
    except ValueError:
        pass
    df      = df_gff_cre.loc[genename]
    mask    = (df[2]=='CDS')
    sub_df = df[mask].reset_index().set_index('transcriptname')
    for ix in set(sub_df.index):
        df_mRNA         = sub_df.loc[ix]
        transcript_name = ix
        
        if isinstance(df_mRNA, pd.Series):
            chromosome = df_mRNA[0]
            left  = df_mRNA[3]
            right = df_mRNA[4]
            echr       = dicChr2N[chromosome]
            covered_array = list(array_contiguity[echr][left-1:right]) # continuity value require minus 1 from right pos
            matched_array = list(match_matrix[echr][left-1:right])
        else:
            chromosome = df_mRNA[0][0]
            array      = df_mRNA[[3,4]].values
            r,c        = np.shape(array)
            if c != 2 :
                print('?!')
                exit()
            covered_array = []
            matched_array = []
            for i in range(r):

                left       = array[i,:][0] #int(df_mRNA[3])
                right      = array[i,:][1] #int(df_mRNA[4])
                echr       = dicChr2N[chromosome]
                contiguity = list(array_contiguity[echr][left-1:right]) # continuity value require minus 1 from right pos
                matched    = list(match_matrix[echr][left-1:right])
                covered_array += contiguity
                matched_array += matched


        covered_array = np.array(covered_array)
        length = len(covered_array)
        dic['mRNA'].append(transcript_name)
        dic['length'].append(length)
        dic['total.depth'].append(sum(covered_array))
        dic['ratio.depth'].append(float(sum(covered_array))/float(length))
        dic['coverage (1x)'].append(len((covered_array >= 1).nonzero()[0])/float(length))
        dic['coverage (10x)'].append(len((covered_array >= 10).nonzero()[0])/float(length))
        dic['coverage (30x)'].append(len((covered_array >= 30).nonzero()[0])/float(length))
        dic['match'].append(sum(matched_array))
        dic['match.ratio'].append(float(sum(matched_array))/float(length))

df_cont = pd.DataFrame(dic)

#mask  = (df_cont['coverage (1x)'] > 0.8) & (df_cont['match.ratio'] > 0.6)
#df_cont_cov = df_cont[mask]
df_cont_cov = df_cont
matrix = df_cont_cov.sort_values(by='total.depth',ascending=False)[['mRNA','coverage (1x)','coverage (10x)','coverage (30x)','ratio.depth','total.depth','match','match.ratio']].values
np.savetxt(file_bam+'.transcripts.all.txt',matrix,fmt='%s',delimiter='\t')
#Outfile = open(file_bam+'.all.txt','w')
#print(file_bam, '\t'.join(map(str,array)),sep='\t',file=Outfile)
