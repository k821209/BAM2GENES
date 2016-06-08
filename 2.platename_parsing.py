#!/usr/bin/python
from __future__ import print_function
filein  = 'all.topthree.txt'
Outfile = open(filein+'.plate.txt','w') 

file_plateinfo = 'plateinfo.txt'

dic = {}
for line in open(file_plateinfo):
    cell = line.strip().split('\t')
    pn   = cell[0]
    wn   = cell[1][0] + '%02d'%int(cell[1][1:])
    dic[(pn,wn)] = cell[2:]
    


for line in open(filein):
    cell = line.strip().split('\t')
    pn = cell[0].split('_')[0].replace('T2-','')
    wn  = cell[0].split('_')[1]  #T2-1_A02_S2.sam.sorted.bam
    print(pn,wn,'\t'.join(cell),'\t'.join(dic[(pn,wn)]),file=Outfile,sep='\t')
