#! /bin/bash

cd ./sources

#extract file
gunzip $1.vcf.gz 

#remove '
sed -i -e s/\'//g $1.vcf
#remove #VAR (because we want to mantain # in headers)
sed -i "s/#VAR/VAR/g" $1.vcf
#remove \x2c which encode for _
sed -i "s/\\\x2c//g" $1.vcf

#sort vcf
vcf-sort $1.vcf > $1_sorted.vcf
#BGzip
bgzip $1_sorted.vcf

#create index .tbi
tabix -p vcf $1_sorted.vcf.gz
#subset chr1 to new file
tabix -h $1_sorted.vcf.gz 1 > $1.vcf
#append remaining chromosomes 
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do tabix $1_sorted.vcf.gz $i >> $1.vcf; done

#remove extra files
rm $1_sorted.vcf.gz
rm $1_sorted.vcf.gz.tbi
#BGzip
bgzip $1.vcf
