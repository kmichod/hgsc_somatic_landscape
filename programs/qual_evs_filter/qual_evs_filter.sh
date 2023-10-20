#!bin/bash
#SBATCH--time=1:00:00
#SBATCH--account=doherty
#SBATCH --partition=redwood
#SBATCH--mail-type=FAIL,END
#SBATCH--mail-user=katherine.lawson-michod@hci.utah.edu

cd /uufs/chpc.utah.edu/common/HIPAA/u1325930/aaces/hgsc_landscape/data/somatic/vcfs/filtered/hg38/moderateconfidence/temp/

for x in $(ls *.vcf) 
do 
vcffilter -f "FILTER = PASS & QUAL > 15" $x > /uufs/chpc.utah.edu/common/HIPAA/u1325930/aaces/hgsc_landscape/data/somatic/vcfs/filtered/hg38/moderateconfidence/${x%%.vcf}.vcf 
done
