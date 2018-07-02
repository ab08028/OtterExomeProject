#! /bin/bash
#$ -cwd
#$ -l h_rt=00:20:00,h_data=1G
#$ -m bea
#$ -N plmxTest
#$ -u ab0808
source /u/local/Modules/default/init/modules.sh
module load python/2.7
plmxDir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/paleomixLinks
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/ancient/paleomixPipeline

$plmxDir/bam_pipeline mkfile > $scriptDir/000_makefile_template.yaml
