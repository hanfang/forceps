#!/usr/bin/env python

### Copyright (c) 2015 Han Fang - hanfang.cshl@gmail.com

### Import necessary libraries
import os
import sys
from scipy import stats
import numpy as np
from numpy.random import randn
import pandas as pd
import datetime
import time


## Require a fasta file as input
if not len(sys.argv) == 3:
    sys.exit("Usage:\tfind_candidates.py <query> <kmer_database> \nError:\trequire two input files, the <query> kmers and the computed <kmer_database>")

### Convert the input to list
Query=sys.argv[1]
Database=sys.argv[2]

pd.set_option('expand_frame_repr', False)

## -----------------
## initial status
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Analysis starts ] \n"


## -----------------
## loading files to start
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Loading files to start ] \n"

query=pd.io.parsers.read_table(Query,header=None,index_col=False)
database=pd.io.parsers.read_table(Database,header=None,index_col=False)
database.columns=['#chr','start','end','qry','strand','ED','MM','INS','DEL','B2B','EDGE','Wobble']

## define dictionaries
kmer={}
min_ed={}

## -----------------
## create db for each kmer
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Creating db for each kmer ] \n"

for i in query.loc[:,0]:
    kmer[i]=pd.DataFrame( database[database.qry==i] )

## -----------------
## find minimum ed for each kmer
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Finding the minimum ed for each kmer ] \n"

for i in query.loc[:,0]:
    min_ed[i]=kmer[i][kmer[i].ED==kmer[i].ED.min()]

## -----------------
## output minimum ed for each kmer
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Export the minimum ED (min ED) for each kmer ] \n"

for i in query.loc[:,0]:
    min_ed[i].to_csv('./min_ed/min_ed.'+i, sep='\t',index=False)

## -----------------
## Concatenating min ED matrix
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Concatenating min ED matrix ] \n"

concat=pd.DataFrame()
for i in query.loc[:,0]:
    concat=concat.append(min_ed[i])

## -----------------
## Calculate the maximum min ED
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Calculating the maximum min ED ] \n"

maxnum_med=concat.ED.max()

## -----------------
## find the kmers with min_ed less than the maximum minimum ed
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Identify the kmers with min ED less than the maximum min ED ] \n"

discarded=concat[concat.ED < maxnum_med].qry

## -----------------
## keep the kmers with max min_ed
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Identify the kmers with maximum min ED ] \n"

max_med=concat[~concat.qry.isin(discarded)]
max_med.to_csv('./min_ed/max_med.txt', sep='\t',index=False)

## -----------------
## remove the kmers with only mm
## -----------------
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print >> sys.stderr, "[",st,"] [ Remove the kmers with only mm ] \n"

onlymm=max_med[max_med.MM == maxnum_med].qry
good_cand=max_med[~max_med.MM.isin(onlymm)]
good_cand.to_csv('./min_ed/good_cand.txt', sep='\t',index=False)

#print 'good_cand'
#print good_cand
    # print concat.loc[:,'qry'] == i
    # if concat.loc[concat.loc[:,'qry'] == i,'ED']==5:
    #print i
    #if concat.loc[concat.loc[:,'qry']==i,:].loc[:,'ED'] <= max_num_med:
    # print concat.loc[concat.loc[:,'qry']==i,:].loc[:,'ED']

'''


cand_list=set(max_med.loc[:,'qry'])
# print cand_list
## exclude the kmers with at mismatches
#for i in max_med.loc[:,'qry']:
    #cand_list.
subset_cand_list=cand_list
for i in cand_list:
    for j in concat.loc[i]
    if concat.loc[i,'MM'] == concat.loc[i,'ED']:
        subset_cand_list.pop(i)

print subset_cand_list
#    no_snps_max_med=max_med.loc[ concat.loc[i:,'MM'] != concat.loc[i:,'ED'] ]

'''
