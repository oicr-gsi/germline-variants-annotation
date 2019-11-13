#!/bin/bash
variants=$1 # file containing path to GATKHaplotypeCaller GVCF files
outputname=$2

logd=`pwd`

script01=`pwd`/01_jointCalling.sh
jointCalledVCF=`pwd`/${outputname}_jointCalled.vcf
jobname01="jointCallVCF_${outputname}"

script02=`pwd`/01_jointGenotyping.sh
jobname02="jointGTVCF_${outputname}"
jointGTVCF=`pwd`/${outputname}_jointGT.vcf

vqsrVCF=`pwd`/${outputname}_jointCalled_recal.vcf
jobname03="vqsr_${outputname}"
script03=`pwd`/03_vqsr.sh

gvcf2vcf=`pwd`/${outputname}_jointCalled_recal.vcf.gz
jobname04="vcfGZ_${outputname}"
script04=`pwd`/04_gvcf2vcf.sh

jobname04b="varmetrics_${outputname}"
script04b=`pwd`/04b_collectVariantMetrics.sh

annotatedVCF=`pwd`/${outputname}_annotated.vep.vcf.gz
jobname05="annotateVCF_${outputname}"
script05=`pwd`/05_annotateVCF.sh


# combined calling
# references
refFasta=/.mounts/labs/PDE/data/reference/hg19/gatk/2.8/ucsc.hg19.fasta
#inputs
#variants
# output
#jointCalledVCF
echo '#!/bin/bash' > $script01
echo "module load java/8" >> $script01
echo "module load gatk/3.6-0" >> $script01
echo "java -jar \${GATK_ROOT}/GenomeAnalysisTK.jar \\" >> $script01
echo "-T CombineGVCFs \\" >> $script01
echo "-R $refFasta \\" >> $script01
for v in `cat $variants`; do
  echo "--variant $v \\" >> $script01
done
echo "-o $jointCalledVCF" >> $script01

# joint GT
# references
refFasta=/.mounts/labs/PDE/data/reference/hg19/gatk/2.8/ucsc.hg19.fasta
#input
#jointCalledVCF
# output
#jointGTVCF
echo '#!/bin/bash' > $script02
echo "module load java/8" >> $script02
echo "module load gatk/3.6-0" >> $script02
echo "java -jar \${GATK_ROOT}/GenomeAnalysisTK.jar \\" >> $script02
echo "-T GenotypeGVCFs \\" >> $script02
echo "-R $refFasta \\" >> $script02
echo "-V $jointCalledVCF \\" >> $script02
echo "-o $jointGTVCF" >> $script02


# VQSR
# references
refFasta=/.mounts/labs/PDE/data/reference/hg19/gatk/2.8/ucsc.hg19.fasta
omniVCF=/.mounts/labs/PDE/data/gatkAnnotationResources/1000G_omni2.5.hg19.sites.vcf
kgSNP=/.mounts/labs/PDE/data/gatkAnnotationResources/1000g_20100804_chr.vcf
hapMapVCF=/.mounts/labs/PDE/data/gatkAnnotationResources/hapmap_3.3.hg19.sites.vcf
dbsnpVCF=/.mounts/labs/PDE/data/gatkAnnotationResources/dbSNP135_chr.vcf
kgIndels=/.mounts/labs/TGL/gsi/databases/1000G_phase1.indels.b37.vcf_chr.vcf
millsIndelsVCF=/.mounts/labs/TGL/gsi/databases/Mills_and_1000G_gold_standard.indels.b37.vcf_chr.vcf
tmp=`pwd`/tmp; mkdir -p $tmp
recalibrate_SNP_recal=$tmp/recalibrate_SNP.recal
recalibrate_SNP_tranches=$tmp/recalibrate_SNP.tranches
recalibrate_SNP_plots_R=$tmp/recalibrate_SNP_plots.R
recalSNPIndelVCF=$tmp/recal.snp.indel.vcf
recalibrate_INDEL_recal=$tmp/recalibrate_INDEL.recal
recalibrate_INDEL_tranches=$tmp/recalibrate_INDEL.tranches
recalibrate_INDEL_plots_R=$tmp/recalibrate_INDEL_plots.R
#input
#jointGTVCF
# output
#vqsrVCF

echo '#!/bin/bash' > $script03
echo "module load java/8" >> $script03
echo "module load gatk/3.6-0" >> $script03
echo "# snp model" >> $script03
echo "java -jar \${GATK_ROOT}/GenomeAnalysisTK.jar \\" >> $script03
echo "    -T VariantRecalibrator \\" >> $script03
echo "    -R $refFasta \\" >> $script03
echo "    -input $jointGTVCF \\" >> $script03
echo "    -U ALLOW_SEQ_DICT_INCOMPATIBILITY \\" >> $script03
echo "    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapMapVCF \\" >> $script03
echo "    --resource:omni,known=false,training=true,truth=true,prior=12.0 $omniVCF \\" >> $script03
echo "    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $kgSNP \\" >> $script03
echo "    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnpVCF \\" >> $script03
echo "    -an DP \\" >> $script03
echo "    -an QD \\" >> $script03
echo "    -an FS \\" >> $script03
echo "    -an SOR \\" >> $script03
echo "    -an MQ \\" >> $script03
echo "    -an MQRankSum \\" >> $script03
echo "    -an ReadPosRankSum \\" >> $script03
echo "    -an InbreedingCoeff \\" >> $script03
echo "    -mode SNP \\" >> $script03
echo "    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\" >> $script03
echo "    -recalFile ${recalibrate_SNP_recal} \\" >> $script03
echo "    -tranchesFile ${recalibrate_SNP_tranches} \\" >> $script03
echo "    -rscriptFile ${recalibrate_SNP_plots_R}" >> $script03
echo "# apply snp model" >> $script03
echo "java -jar \${GATK_ROOT}/GenomeAnalysisTK.jar \\" >> $script03
echo "    -T ApplyRecalibration \\" >> $script03
echo "    -R /.mounts/labs/PDE/data/reference/hg19/gatk/2.8/ucsc.hg19.fasta \\" >> $script03
echo "    -input $jointGTVCF \\" >> $script03
echo "    -mode SNP \\" >> $script03
echo "    --ts_filter_level 99.0 \\" >> $script03
echo "    -recalFile ${recalibrate_SNP_recal} \\" >> $script03
echo "    -tranchesFile ${recalibrate_SNP_tranches} \\" >> $script03
echo "    -o $recalSNPIndelVCF" >> $script03
echo "#build indel model" >> $script03
echo "java -jar \${GATK_ROOT}/GenomeAnalysisTK.jar \\" >> $script03
echo "    -T VariantRecalibrator \\" >> $script03
echo "    -R /.mounts/labs/PDE/data/reference/hg19/gatk/2.8/ucsc.hg19.fasta \\" >> $script03
echo "    -input $recalSNPIndelVCF \\" >> $script03
echo "    --resource:mills,known=false,training=true,truth=true,prior=12.0 $millsIndelsVCF \\" >> $script03
echo "    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $kgIndels \\" >> $script03
echo "    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnpVCF \\" >> $script03
echo "    -U ALLOW_SEQ_DICT_INCOMPATIBILITY \\" >> $script03
echo "    -an QD \\" >> $script03
echo "    -an DP \\" >> $script03
echo "    -an FS \\" >> $script03
echo "    -an SOR \\" >> $script03
echo "    -an MQRankSum \\" >> $script03
echo "    -an ReadPosRankSum \\" >> $script03
echo "    -an InbreedingCoeff \\" >> $script03
echo "    -mode INDEL \\" >> $script03
echo "    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\" >> $script03
echo "    --maxGaussians 4 \\" >> $script03
echo "    -recalFile ${recalibrate_INDEL_recal} \\" >> $script03
echo "    -tranchesFile ${recalibrate_INDEL_tranches} \\" >> $script03
echo "    -rscriptFile ${recalibrate_INDEL_plots_R} " >> $script03
echo "# apply indel model " >> $script03
echo "java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar \\" >> $script03
echo "    -T ApplyRecalibration \\" >> $script03
echo "    -R $refFasta \\" >> $script03
echo "    -input $recalSNPIndelVCF \\" >> $script03
echo "    -mode INDEL \\" >> $script03
echo "    --ts_filter_level 99.0 \\" >> $script03
echo "    -recalFile ${recalibrate_INDEL_recal} \\" >> $script03
echo "    -tranchesFile ${recalibrate_INDEL_tranches} \\" >> $script03
echo "    -o $vqsrVCF " >> $script03

# gvcf2vcf
#input
#vqsrVCF
# output
#gvcf2vcf
echo '#!/bin/bash' > $script04
echo "module load tabix/1.9; cat $vqsrVCF | grep -v END | grep -v GVCFBlock | sed 's/,<NON_REF>//g' | bgzip -c > $gvcf2vcf; tabix -p vcf $gvcf2vcf" >> $script04

# run some variant stats
echo '#!/bin/bash' > $script04b
echo "module load picard/2.19.2 " >> $script04b
echo "ava -jar ${PICARD_ROOT}/picard.jar CollectVariantCallingMetrics \\" >> $script04b
echo "INPUT=$vqsrVCF \\" >> $script04b
echo "OUTPUT=${vqsrVCF}.metrics \\" >> $script04b
echo "DBSNP=$dbsnpVCF " >> $script04b


# VEP annotate
# references
DBNSFP=/.mounts/labs/TGL/gsi/databases/vep_plugin_data/dbNSFP_hg19.gz
FASTA=/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19/Genomic/references/fasta/hg19.fa
LOFSCORE=/.mounts/labs/TGL/gsi/databases/vep_plugin_data/LoFtool_scores.txt
#input
#gvcf2vcf
# output
#annotatedVCF
echo '#!/bin/bash' > $script05
echo "DBNSFP=$DBNSFP" >> $script05
echo "FASTA=$FASTA" >> $script05
echo "LOFSCORE=$LOFSCORE" >> $script05
echo "# perl/5.22.2-tgl" >> $script05
echo "export LD_LIBRARY_PATH=/oicr/local/analysis/sw/perl/perl-5.22.2-tgl/lib:$LD_LIBRARY_PATH;" >> $script05
echo "export PERL5LIB=/oicr/local/analysis/sw/perl/perl-5.22.2-tgl/lib:$PERL5LIB" >> $script05
echo "export PATH=/oicr/local/analysis/sw/perl/perl-5.22.2-tgl/bin:$PATH" >> $script05
echo "# vep/92" >> $script05
echo "export PATH=/oicr/local/analysis/sw/vep/vep92:$PATH" >> $script05
echo "export PATH=/oicr/local/analysis/sw/vep/vep92/htslib:$PATH" >> $script05
echo "export PATH=/oicr/local/analysis/sw/vep/vep92/samtools/bin:$PATH" >> $script05
echo "export PERL5LIB=/oicr/local/analysis/sw/vep/vep92:$PATH" >> $script05
echo "export VEP_PATH=/oicr/local/analysis/sw/vep/vep92;" >> $script05
echo "export VEP_DATA=/oicr/local/analysis/sw/vep/vep92/.cache;" >> $script05
echo "# vcf2maf" >> $script05
echo "export PATH=/.mounts/labs/PDE/Modules/sw/vcf2maf/mskcc-vcf2maf-decbf60:$PATH;" >> $script05
echo "vep=/oicr/local/analysis/sw/vep/vep92/vep" >> $script05
echo "$vep -i $gvcf2vcf \\" >> $script05
echo "                --offline --dir_cache $VEP_DATA \\" >> $script05
echo "                --assembly GRCh37 \\" >> $script05
echo "                --database $VEP_DATA \\" >> $script05
echo "                --force_overwrite \\" >> $script05
echo "                --fasta $FASTA \\" >> $script05
echo "                --variant_class \\" >> $script05
echo "                --hgvs \\" >> $script05
echo "                --symbol \\" >> $script05
echo "                --canonical \\" >> $script05
echo "                --check_existing \\" >> $script05
echo "                --humdiv \\" >> $script05
echo "                --sift b \\" >> $script05
echo "                --polyphen b \\" >> $script05
echo "                --plugin dbNSFP,$DBNSFP,genename,clinvar_golden_stars,clinvar_clnsig,1000Gp3_AF,ExAC_AF,gnomAD_exomes_AF,gnomAD_genomes_AF,SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred,FATHMM_pred,REVEL_score,CADD_phred,GERP++_RS \\" >> $script05
echo "                --plugin LoFtool,$LOFSCORE \\" >> $script05
echo "                --vcf -o $annotatedVCF" >> $script05

# clean up

chmod +x $script01
qsub -V -l h_vmem=32g -N $jobname01 -e $logd -o $logd $script01
chmod +x $script02
qsub -V -l h_vmem=32g -hold_jid $jobname01 -N $jobname02 -e $logd -o $logd $script02
chmod +x $script03
qsub -V -l h_vmem=32g -hold_jid $jobname02 -N $jobname03 -e $logd -o $logd $script03
chmod +x $script04
qsub -V -l h_vmem=32g -hold_jid $jobname03 -N $jobname04 -e $logd -o $logd $script04
chmod +x $script04b
qsub -V -l h_vmem=32g -hold_jid $jobname03 -N $jobname04b -e $logd -o $logd $script04b
chmod +x $script05
qsub -V -l h_vmem=32g -hold_jid $jobname04 -N $jobname05 -e $logd -o $logd $script05
