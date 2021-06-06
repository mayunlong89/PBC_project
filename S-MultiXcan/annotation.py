
#!usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 15:57:51 2021

@author: Yunlong Ma
"""

import os

os.getcwd()

#f1 = open("file1.txt","r+")
#f2 = open("file2.txt","r+")
#f3 = open("results2.txt","w+")

f1 = open("PBC_GWAS_UKBiobank_summary","r+")
f2 = open("PBC_SNP_CHR_POS.txt","r+")
f3 = open("PBC_GWAS_UKBiobank_summary_final","w+")


data = []
for line in f1:
    dd = line.strip().split()
    data.append(dd)


anno_file =[]
for line2 in f2:
    aa = line2.strip().split()
    anno_file.append(aa)
#print(anno_file)
#print(data)


for i in range(0,len(data)):
    data_anno =[]
    for j in range(0,len(anno_file)):
        if data[i][0] == anno_file[j][0]:
            CHR = anno_file[j][1]
            POS = anno_file[j][2]
            DD = "\t".join(data[i])+'\t'+ CHR + '\t' + POS +"\n"
    data_anno.append(DD)
    #anno = "\t".join(data[i])+'\t'+ DD+"\n"
    anno = "\t".join(data_anno)
    f3.write(anno)


