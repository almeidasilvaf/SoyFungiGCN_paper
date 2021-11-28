#!/usr/bin/bash

cd ../products/result_files
[ ! -f soysnp50k.vcf.gz ] && curl -k -o soysnp50k.vcf.gz https://soybase.org/snps/soysnp50k_wm82.a2_41317.vcf.gz
[ ! -f soysnp50k.vcf ] && gunzip soysnp50k.vcf.gz
cut -f1 snps_and_effects.txt | grep -v SNP > ids.txt

# Extract indices of rows that match
touch filtered_soysnp50k.txt
for i in $(cat ids.txt)
	do grep "$i" soysnp50k.vcf >> filtered_soysnp50k.txt
done
cat soysnp50k.vcf | grep PI > header.txt
cat header.txt filtered_soysnp50k.txt > filtered_soysnp50k_final.txt

rm soysnp50k.vcf
