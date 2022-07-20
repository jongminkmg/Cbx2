#### Purpose:
#### To make smaller fastq file with only R1 and R2 sequences
#### for faster processing.
####
#### Input: R1-R2 combined 8-line repeating fastq
#### e.g.) head S07-08_R12.txt
#### @VH00377:2:AAAGMCVHV:1:1101:44457:1019 1:N:0:GTAGCCCTGT+ATAGATGCTC
#### @VH00377:2:AAAGMCVHV:1:1101:44457:1019 2:N:0:GTAGCCCTGT+ATAGATGCTC
#### CTCGGAGTCATCATTCTCTTGGTTGG
#### TAGAGGACGAACTGCTGGATTTGGATTTGGATGGCGCATCTGGTTCTTGGAGGACCAGCCGCGCCACTTGACCAGGTACTCCAGCTTGCC
#### +
#### +
#### CCCCCCCCCCCCCCCCCCCCCCCCCC
#### CCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;;CCC
####
#### Output: Return only line3(R1 sequence), line4(R2 sequence)
####
#### By Jongmin Kim (jongminkmg@gmail.com, 2021. 09. 26)



File2=open('S07-08_trimmed.txt', mode='w')

import csv
import codecs

counter = 1

with open('S07-08_R12.txt', 'rb') as f1:
        reader_Unique = csv.reader(codecs.iterdecode(f1, 'utf-8'), delimiter='\t')
        for row_Unique in reader_Unique:
            if counter == 3:
                File2.write(str(row_Unique)[2:-2]+'\n')
            if counter == 4:
                File2.write(str(row_Unique)[2:-2]+'\n')
            counter = counter+ 1
            
            if counter == 9:
                counter = 1

File2.close()
