#!bin/bash

cd /hgsc_somatic_landscape/data/somatic/vcfs/filtered/hg38/highconfidence/temp/

for x in $(ls *.vcf) 
do 
vcffilter -f "FILTER = PASS & QUAL > 15" $x > hgsc_somatic_landscape/data/somatic/vcfs/filtered/hg38/highconfidence/${x%%.vcf}.vcf 
done
