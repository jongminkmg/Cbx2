#### Purpose:
#### Obtain the number of Cbx2 WT 'unique(UMI based)' & Cbx2 mutant (tm1d) reads
#### from Cbx2 amplicon reads (5' single cell 10X read structure). 

#### For input & executed output files:

#### Input1: text file with barcode sequences.
#### $ head S01_Cluster2_KMeans7.txt
#### AAACGGGAGCCACTAT
#### AAACGGGTCCAATGGT
#### AAAGATGAGCGTTCCG

#### Input2: trimmed fastq file
#### R1 \t R2
#### $ head S05-06_trimmed3half.txt
#### GGTGAAGGTCGATTGTTCAGCGGGGC      ACAGAGGACGAACTGCTGGATTTGGATTTGGATGGCGCATCTGGTTCCTTGAGCTTGGAGCGCCGGCTGCAGGAGGATGTGACTGTGTGT
#### CTCGGAGTCATCATTCTCTTGGTTGG      TAGAGGACGAACTGCTGGATTTGGATTTGGATGGCGCATCTGGTTCTTGGAGGACCAGCCGCGCCACTTGACCAGGTACTCCAGCTTGCC
#### AGTGTCAGTGTTTGGTAGCACCGTGC      TAGAGGACGAACTGCTGGATTTGGATTTGGATGGCGCATCTGGTTCTTGGAGGACCAGCCGCGCCACTTGACCAGGTACTCCAGCTTGCC

#### Output1: 
#### Line1: All UMI sequences
#### Line2: Number of WT called
#### Line3: Number of Mut called

#### Output2: 
#### Barcode \t Number of WT \t Number of Mutant

#### Last updated: 2021-09-28, Jongmin Kim (jongminkmg@gmail.com)




File2=open('PO_Seurat_intAB_S03_Comprehensive.txt', mode='w')
File3=open('PO_Seurat_intAB_S03_OnlyCounts.txt', mode='w')

import csv
import codecs

counter = 1
counter2 = 0

barcode = 'null'

genotype = -9 # default -9

CountWT = 0
CountMut = 0

R1 =[]
WT = []
Mut = []



with open('Seurat_intAB_S03.txt', 'rb') as f1:
    reader_barcode = csv.reader(codecs.iterdecode(f1, 'utf-8'), delimiter='\t')
    for row_barcode in reader_barcode:
        barcode = str(row_barcode)[2:-2]
             
        with open('S07-08_trimmed2.txt', 'rb') as f2:    
            reader_Unique = csv.reader(codecs.iterdecode(f2, 'utf-8'), delimiter='\t')
            for row_Unique in reader_Unique:            
                if row_Unique[0].find(barcode) > -1:

                    try:
                        SavedIndex = R1.index(row_Unique[0])
                        if row_Unique[1].find('TGGCGCATCTGGTTCCTTGAGCTTGGAGCG') > -1:
                            WT[SavedIndex] += 1
                        elif row_Unique[1].find('TGGCGCATCTGGTTCTTGGAGGACCAGCCG') > -1:
                            Mut[SavedIndex] += 1
                        
                    except ValueError:
                        R1.append(row_Unique[0])
                        NewIndex = len(R1)-1
                        if row_Unique[1].find('TGGCGCATCTGGTTCCTTGAGCTTGGAGCG') > -1:
                            WT.append(1)
                            Mut.append(0)
                        elif row_Unique[1].find('TGGCGCATCTGGTTCTTGGAGGACCAGCCG') > -1:
                            WT.append(0)
                            Mut.append(1)
                        else:
                            WT.append(0)
                            Mut.append(0)

            print ('WT: ', WT, '\n')
            print ('Mut: ', Mut, '\n')
            for i in range(len(R1)):
                File2.write(R1[i]+'\t')
            File2.write('\n')
            for i in range(len(WT)):
                File2.write(str(WT[i])+'\t')
            File2.write('\n')
            for i in range(len(Mut)):
                File2.write(str(Mut[i])+'\t')
            File2.write('\n\n')

            for i in range(len(R1)):
                if WT[i] == 0 and Mut[i] > 0:
                    CountMut += 1
                elif WT[i] > 0 and Mut[i] == 0:
                    CountWT += 1
                elif WT[i] >0 and Mut[i] == 0:
                    if WT[i]/Mut[i] > 10:
                        CountWT += 1
                    elif Mut[i]/WT[i] >10:
                        CountMut += 1
                        
            print ('row_barcode: ', barcode, CountWT, CountMut)
            File3.write(barcode+'\t'+str(CountWT)+'\t'+str(CountMut)+'\n')
            R1 =[]
            WT = []
            Mut = []
            CountWT = 0
            CountMut = 0


File2.close()
File3.close()
