#!/bin/bash

## this file runs the VCF annotator and parser which will filter called somatic variants based on specific criteria

## read depth -d, minimum tumor AF -m, minimum alternative observations -w, and BKZ score -z set based on criteria for low, medium, and high confidence criteria

java -jar ~/USeq/Apps/AnnotatedVcfParser -v ~/hgsc_somatic_landscape/data -s ~/hgsc_somatic_landscape/data/somatic_highconfidence -d 30 -m 0.15 -w 4 -z 3

## after vcf parser remove any VCFs failing these parameters as well as the hyper-mutated sample 

cd ~/hgsc_somatic_landscape/data/somatic_highconfidence

rm *Fail*

## Unzip VCF files passing VCF parser

gunzip *.vcf.gz

## Rename files 

for filename in *vcf ; 
do replacement=$(echo $filename | sed -e 's,^.*/,,' -e 's,\.[^\.]*$,,')
[ -f "$filename" ] && sed -i "s,TUMOR,$replacement," "$filename"
done

## and trim file names 
for x in *vcf; do 
echo $x; grep -rl '_Anno_Hg38.anno_Pass' $x | xargs sed -i 's/_Anno_Hg38.anno_Pass//g'; 
done

