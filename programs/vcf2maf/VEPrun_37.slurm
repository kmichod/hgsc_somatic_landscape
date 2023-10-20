#!/bin/bash

#SBATCH --account=doherty

#SBATCH --partition=redwood

#SBATCH -N 1

#SBATCH -t 72:00:00

#Lawson-Michod


set -e

module load vep samtools

start=$(date +'%s'); rm -f FAILED COMPLETE QUEUED; touch STARTED

echo -e "\n---------- Starting -------- $((($(date +'%s') - $start)/60)) min"



for x in $(ls *Pass_moderate.vcf)

do

echo $x

/uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/Modules/VEP106/run_37maf.sh --ref-fasta /uufs/chpc.utah.edu/common/HIPAA/u1325930/aaces/hgsc_landscape/reference/hg19/hg19_ucsc.hg19.fasta --input-vcf $x --ncbi-build GRCh37 --retain-info T_DP,T_AF,N_DP,N_AF,BKZ,VCFSS,QSI,SomaticEVS,QSS,SNVSB,ANN,LOF,NMD,DBVARID,ALLELEID,CLNSIG,CLNVCSO,OLD_VARIANT,CLNREVSTAT,RS,CLNDNINCL,ORIGIN,MC,CLNDN,CLNVC,CLNVI,AF_EXAC,OLD_CLUMPED,AF_ESP,CLNSIGINCL,CLNDISDB,GENEINFO,CLNDISDBINCL,AF_TGP,CLNSIGCONF,SSR,OLD_MULTIALLELIC --tumor-id ${x%%_Illumina*} --normal-id Normal --output-maf ${x%%.vcf}.maf --maf-center HCI --vep-forks 4

#--remap-chain /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/Modules/VEP106/data/hg19_to_GRCh37.chain

done



echo -e "\n---------- Complete! -------- $((($(date +'%s') - $start)/60)) min total"



touch COMPLETE
