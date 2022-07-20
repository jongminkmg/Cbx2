#### Purpose:
#### 1) Remove R1-R2 pair when R2 sequence is GGGGG...GGGGG
#### 2) Rearrange sequence files so that it will be R1,'\t',R2
####
#### Input: R1(line1)-R2(line2) only sequence file
#### e.g.)
#### CTCGGAGTCATCATTCTCTTGGTTGG
#### TAGAGGACGAACTGCTGGATTTGGATTTGGATGGCGCATCTGGTTCTTGGAGGACCAGCCGCGCCACTTGACCAGGTACTCCAGCTTGC
#### .
#### .
#### Output: R1,'\t',R2,'\n' (w/o 'GGGGGGGGGGG' wierd sequences)
####
#### Jongmin Kim (jongminkmg@gmail.com, 2021. 09. 26)



File2=open('S07-08_trimmed2.txt', mode='w')

import csv
import codecs

counter = 1
R1 = 'NA'
R2 = 'NA'

with open('S07-08_trimmed.txt', 'rb') as f1:
        reader_Unique = csv.reader(codecs.iterdecode(f1, 'utf-8'), delimiter='\t')
        for row_Unique in reader_Unique:
            if counter == 1:
                R1 = str(row_Unique)[2:-2]
                
            if counter == 2:
                R2 = str(row_Unique)[2:-2]
                if R2.find('GGGGGGGGGGGGGGGGGGGG') == -1:
                    File2.write(R1+'\t'+R2+'\n')
                
            counter = counter+ 1
            
            if counter == 3:
                counter = 1
                R1 = 'NA'
                R2 = 'NA'

File2.close()
