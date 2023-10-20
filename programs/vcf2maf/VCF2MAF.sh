#!/bin/bash

set -e

module load vep samtools

start=$(date +'%s'); rm -f FAILED COMPLETE QUEUED; touch STARTED

echo -e "\n---------- Starting -------- $((($(date +'%s') - $start)/60)) min"

for x in $(ls *.vcf)

do

echo $x

/uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/Modules/VEP106/run_38maf.sh \ --ref-fasta /uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/TNRunner/GATKResourceBundleAug2021/Homo_sapiens_assembly38.fasta \ --input-vcf $x --ncbi-build GRCh38 --retain-info T_DP,T_AF,N_DP,N_AF,BKZ,VCFSS,QSI,SomaticEVS,QSS,SNVSB,ANN,LOF,NMD,DBVARID,ALLELEID,CLNSIG,CLNVCSO,OLD_VARIANT,CLNREVSTAT,RS,CLNDNINCL,ORIGIN,MC,CLNDN,CLNVC,CLNVI,AF_EXAC,OLD_CLUMPED,AF_ESP,CLNSIGINCL,CLNDISDB,GENEINFO,CLNDISDBINCL,AF_TGP,CLNSIGCONF,SSR,OLD_MULTIALLELIC \ --tumor-id ${x%%_Illumina*} --normal-id Normal --output-maf ${x%%.vcf}.maf \ --maf-center HCI --vep-forks 4 \

#--remap-chain /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/Modules/VEP106/data/hg38_to_GRCh38.chain

done



echo -e "\n---------- Complete! -------- $((($(date +'%s') - $start)/60)) min total"



touch COMPLETE
