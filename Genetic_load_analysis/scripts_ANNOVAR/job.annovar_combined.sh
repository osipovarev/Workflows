#!/bin/bash
#SBATCH -J WGR_ANNOVAR
#SBATCH -o ./logs/ANNOVAR/WGR_ANNOVAR_%a_%A.log
#SBATCH -e ./logs/ANNOVAR/WGR_ANNOVAR_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 24:00:00
#SBATCH -c 32  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=32000  # Requested Memory
#SBATCH --mail-type=END


WDIR=/nese/meclab/Katya/ANNOVAR_turtles/
DB=/nese/meclab/Katya/Scripts/annovar/turtledb/
VCF=/nese/meclab/Lisa/Dc_WGR_analyses_docsMay2024/Final_filtered_combvcf/DerCor_combined_filtered.vcf
BUILD=derCor
PROTOCOL=genes

table_annovar.pl $VCF $DB -out derCor_annovar -buildver $BUILD -protocol $PROTOCOL -operation g -remove -polish -vcfinput -nastring .

