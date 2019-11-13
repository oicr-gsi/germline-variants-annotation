#!/bin/bash
script01=`pwd`/01_jointCalling.sh
logd=`pwd`

variants=$1 # file containing path to GATKHaplotypeCaller GVCF files
outputname=$2
jointCalledVCF=`pwd`/${outputname}_jointCalled.vcf
jobname01="jointCallVCF"

script02=`pwd`/01_jointGenotyping.sh
jobname02="jointGTVCF"
jointGTVCF=`pwd`/${outputname}_jointGT.vcf

gvcf2vcf=`pwd`/${outputname}_jointGT.vcf.gz
jobname03="vcfGZ"
script03=`pwd`/03_gvcf2vcf.sh

annotatedVCF=`pwd`/${outputname}annotated.vep.vcf.gz
jobname04="annotateVCF"
script04=`pwd`/04_annotateVCF.sh


echo "#!/bin/bash" > $script01
echo "module load java/8" >> $script01
echo "module load gatk/3.6-0" >> $script01

echo "java -jar /.mounts/labs/PDE/modulator/sw/Debian8.11/gatk-3.6-0/GenomeAnalysisTK.jar \\" >> $script01
echo "-T CombineGVCFs \\" >> $script01
echo "-R /.mounts/labs/PDE/data/reference/hg19/gatk/2.8/ucsc.hg19.fasta \\" >> $script01
for v in `cat $variants`; do
  echo "--variant $v \\" >> $script01
done
echo "-o $jointCalledVCF" >> $script01



echo "#!/bin/bash" > $script02
echo "module load java/8" >> $script02
echo "module load gatk/3.6-0" >> $script02
echo "java -jar /.mounts/labs/PDE/modulator/sw/Debian8.11/gatk-3.6-0/GenomeAnalysisTK.jar \\" >> $script02
echo "-T GenotypeGVCFs \\" >> $script02
echo "-R /.mounts/labs/PDE/data/reference/hg19/gatk/2.8/ucsc.hg19.fasta \\" >> $script02
echo "-V $jointCalledVCF \\">> $script02
echo "-o $jointGTVCF" >> $script02



echo "#!/bin/bash" > $script03
echo "module load tabix/1.9; cat $jointGTVCF | grep -v END | grep -v GVCFBlock | sed 's/,<NON_REF>//g' | bgzip -c > $gvcf2vcf; tabix -p vcf $gvcf2vcf" >> $script03



echo "#!/bin/bash" > $script04
echo "DBNSFP=/.mounts/labs/TGL/gsi/databases/vep_plugin_data/dbNSFP_hg19.gz" >> $script04
echo "FASTA=/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19/Genomic/references/fasta/hg19.fa" >> $script04
echo "LOFSCORE=/.mounts/labs/TGL/gsi/databases/vep_plugin_data/LoFtool_scores.txt" >> $script04
echo "# perl/5.22.2-tgl" >> $script04
echo "export LD_LIBRARY_PATH=/oicr/local/analysis/sw/perl/perl-5.22.2-tgl/lib:$LD_LIBRARY_PATH;" >> $script04
echo "export PERL5LIB=/oicr/local/analysis/sw/perl/perl-5.22.2-tgl/lib:$PERL5LIB" >> $script04
echo "export PATH=/oicr/local/analysis/sw/perl/perl-5.22.2-tgl/bin:$PATH" >> $script04
echo "# vep/92" >> $script04
echo "export PATH=/oicr/local/analysis/sw/vep/vep92:$PATH" >> $script04
echo "export PATH=/oicr/local/analysis/sw/vep/vep92/htslib:$PATH" >> $script04
echo "export PATH=/oicr/local/analysis/sw/vep/vep92/samtools/bin:$PATH" >> $script04
echo "export PERL5LIB=/oicr/local/analysis/sw/vep/vep92:$PATH" >> $script04
echo "export VEP_PATH=/oicr/local/analysis/sw/vep/vep92;" >> $script04
echo "export VEP_DATA=/oicr/local/analysis/sw/vep/vep92/.cache;" >> $script04
echo "# vcf2maf" >> $script04
echo "export PATH=/.mounts/labs/PDE/Modules/sw/vcf2maf/mskcc-vcf2maf-decbf60:$PATH;" >> $script04
# echo "CURRDIR=/.mounts/labs/TGL/gsi/pipeline/data/TGL22/EXOME/Seqware_GATKHaplotypeCaller/JointCalling/
echo "vep=/oicr/local/analysis/sw/vep/vep92/vep" >> $script04
echo "$vep -i $gvcf2vcf \\" >> $script04
echo "                --offline --dir_cache $VEP_DATA \\" >> $script04
echo "                --assembly GRCh37 \\" >> $script04
echo "                --database $VEP_DATA \\" >> $script04
echo "                --force_overwrite \\" >> $script04
echo "                --fasta $FASTA \\" >> $script04
echo "                --variant_class \\" >> $script04
echo "                --hgvs \\" >> $script04
echo "                --symbol \\" >> $script04
echo "                --canonical \\" >> $script04
echo "                --check_existing \\" >> $script04
echo "                --humdiv \\" >> $script04
echo "                --sift b \\" >> $script04
echo "                --polyphen b \\" >> $script04
echo "                --plugin dbNSFP,/.mounts/labs/TGL/gsi/databases/vep_plugin_data/dbNSFP_hg19.gz,genename,clinvar_golden_stars,clinvar_clnsig,1000Gp3_AF,ExAC_AF,gnomAD_exomes_AF,gnomAD_genomes_AF,SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred,FATHMM_pred,REVEL_score,CADD_phred,GERP++_RS \\" >> $script04
echo "                --plugin LoFtool,$LOFSCORE \\" >> $script04
echo "                --vcf -o $annotatedVCF" >> $script04


chmod +x $script01
qsub -V -l h_vmem=32g -N $jobname01 -e $logd -o $logd $script01
chmod +x $script02
qsub -V -l h_vmem=32g -hold_jid $jobname01 -N $jobname02 -e $logd -o $logd $script02
chmod +x $script03
qsub -V -l h_vmem=32g -hold_jid $jobname02 -N $jobname03 -e $logd -o $logd $script03
chmod +x $script04
qsub -V -l h_vmem=32g -hold_jid $jobname03 -N $jobname04 -e $logd -o $logd $script04
