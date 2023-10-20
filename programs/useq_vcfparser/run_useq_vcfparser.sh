#!/uqsr/bin/bash


## this file runs the VCF annotator and parser which will filter called somatic variants based on specific criteria

## read depth -d, minimum tumor AF -m, minimum alternative observations -w, and BKZ score -z set based on criteria for low, medium, and high confidence criteria

RD = ##input minimum read depth
MAF = ##input minimum allele frequency
ALT = ##input minimum alt observations
BKZ = ##input minimum BKZ score 

java -jar ./USeq/Apps/AnnotatedVcfParser -v ./data/somatic/vcfs/master/* -s ./data/somatic/vcfs/filtered/hg38/highconfidence/temp/* -d $RD -m $MAF -w $ALT -z $BKZ

## after vcf parser remove any VCFs failing these parameters as well as the hyper-mutated sample 

cd ./data/somatic/vcfs/filtered/hg38/highconfidence/temp/

rm *Fail*
rm *197-150015*

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

