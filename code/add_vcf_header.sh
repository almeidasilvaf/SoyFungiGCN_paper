
###################
#
# Code to download soybean .VCF file and create a file with VCF headers
#
###################

cd ../products/result_files
[ ! -f soysnp50k.vcf.gz ] && curl -k -o soysnp50k.vcf.gz https://soybase.org/snps/soysnp50k_wm82.a2_41317.vcf.gz
