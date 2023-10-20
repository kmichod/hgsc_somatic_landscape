#!/bin/bash
#SBATCH--time=10:00:00
#SBATCH--nodes=1
#SBATCH--ntasks=1
#SBATCH--account=doherty
#SBATCH --partition=redwood
#SBATCH--mail-type=FAIL,END
#SBATCH--mail-user=katherine.lawson-michod@hci.utah.edu

## this file runs the VCF annotator and parser which will filter called somatic variants based on specific criteria

## read depth -d, minimum tumor AF -m, minimum alternative observations -w, and BKZ score -z set based on criteria for low, medium, and high confidence criteria

java -jar /uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/TNRunner/BioApps/USeq/Apps/AnnotatedVcfParser -v /uufs/chpc.utah.edu/common/HIPAA/u1325930/aaces/hgsc_landscape/data/somatic/vcfs/sb_master/ -s /uufs/chpc.utah.edu/common/HIPAA/u1325930/aaces/hgsc_landscape/data/somatic/vcfs/filtered/hg38/moderateconfidence/temp/ -d 25 -m .075 -w 4 -z 3

## after vcf parser remove any VCFs failing these parameters as well as the hyper-mutated sample 

cd /uufs/chpc.utah.edu/common/HIPAA/u1325930/aaces/hgsc_landscape/data/somatic/vcfs/filtered/hg38/moderateconfidence/temp/

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

