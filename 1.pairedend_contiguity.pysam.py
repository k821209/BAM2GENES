from __future__ import print_function
import subprocess
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
sys.path.append('/ref/analysis/pipelines/')
import pysam
import kang
file_bam = 'intron3000.merge.sorted.bam'
file_fa = 'Creinhardtii_281_v5.0.fa'
dicHD2seq = kang.Fasta2dic(file_fa)

#chromosome, left, right = 'chromosome_1', 1000, 20000

rows              = len(dicHD2seq.keys())
chromosomes       = dicHD2seq.keys()
chromosomes.sort()
dicN2chr          = dict(enumerate(chromosomes))
dicChr2N          = {b:a for a,b in dicN2chr.iteritems()}
columns           = max([len(x) for x in dicHD2seq.values()])-1
continuity_matrix = np.zeros([rows,columns],dtype=np.int)
Outfile = open('chromosome.map.txt','w')
for a,b in dicChr2N.iteritems():
    print(a,b,sep='\t',file=Outfile)


print('start loop')

samfile = pysam.Samfile( file_bam, "rb" )
it      = samfile.fetch()
for line in tqdm(it):#$open('temp.sam.cut'): # should be changed to zero base map
    # Check qual
    if line.is_proper_pair == False:
        continue
    if line.is_duplicate   == True:
        continue
    if line.is_qcfail      == True:
        continue
    if line.is_secondary   == True:
        continue
    # Check qual end
    
    chromosome   = line.reference_name
    startpos     = line.reference_start # zero based
    qname        = line.qname
    echr         = dicChr2N[chromosome]
    fragmentsize = line.tlen
    if line.mpos - startpos > 0 : 
        continuity_matrix[echr,startpos:line.mpos] += 1  # list characteristic can utillize fragment size itself.
    else:
        continuity_matrix[echr,startpos:startpos+line.reference_length-1] += 1

print('saving..')
np.save('%s.paired_end_continuity.np'%file_bam,continuity_matrix)

#np.savetxt('test.txt',continuity_matrix)    


