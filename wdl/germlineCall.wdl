version 1.0

workflow  HaplotypeCallerVariantCalling {

  input {
    File varFile
    String outputFileNamePrefix = outputPrefix
    # references --> to be replaced by data modulator
    String? refFasta
    String? omniVCF
    String? kgSNP
    String? hapMapVCF
    String? dbsnpVCF
    String? kgIndels
    String? millsIndelsVCF
    String? DBNSFP
    String? FASTA
    String? LOFSCORE
    String? tmp # temp directory
  }

  call combineGVCF {
    input:
      vcfs=variantfile,
      outputFileNamePrefix=outputFileNamePrefix,
      refFasta=refFasta
 }

  call jointGenotype {
    input:
      outputFileNamePrefix=outputFileNamePrefix,
      refFasta=refFasta,
      combinedGCVF=combineGVCF.combinedGVCF

  }

  call vqsr {
    input:
      outputFileNamePrefix=outputFileNamePrefix,
      refFasta=refFasta,
      jointGenotypedVCF=jointGenotype.jointGenotypedVCF,
      tmp=tmp

  }

  call varmetrics {
    input:
      outputFileNamePrefix=outputFileNamePrefix,
      refFasta=refFasta,
      vqsrVCF=vqsr.vqsrVCF

  }

  call annotateVCF {
    input:
      outputFileNamePrefix=outputFileNamePrefix,
      vqsrVCF=vqsr.vqsrVCF,
      omniVCF=omniVCF,
      kgSNP=kgSNP,
      hapMapVCF=hapMapVCF,
      dbsnpVCF=dbsnpVCF,
      kgIndels=kgIndels,
      millsIndelsVCF=millsIndelsVCF,
      DBNSFP=DBNSFP,
      FASTA=FASTA,
      LOFSCORE=LOFSCORE
  }

  output {
    File varMetrics = varmetrics.varMetrics
    File annotatedVCF = annotateVCF.annotatedVCF
  }

  parameter_meta {
    variantsfile: "Input bam."
    bamIndex: "Input bam index (must be .bam.bai)."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    windowSize: "The size of non-overlapping windows."
    minimumMappingQuality: "Mapping quality value below which reads are ignored."
    chromosomesToAnalyze: "Chromosomes in the bam reference file."
  }

  meta {
    author: "Prisni Rath"
    email: "Prisni.Rath@oicr.on.ca"
    description: "Workflow for combining GVCFs from haplotype caller"
    dependencies: [
      {
        name: "java/8",
        url: "https://github.com/AdoptOpenJDK/openjdk8-upstream-binaries/releases/download/jdk8u222-b10/OpenJDK8U-jdk_x64_linux_8u222b10.tar.gz"
      },
      {
        name: "gatk/3.6-0",
        url: "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.6-0-g89b7209"
      }
    ]
  }

}

task combineGVCF {
  input{
    File vcfs,
    String outputFileNamePrefix,
    File refFasta,
    String? modules = "java/8,gatk/3.6-0",
    Int? mem = 32
  }

  command <<<
    set -o pipefail

    vars=`cat ~{vcfs}`
    varString=`for v in $vars; do echo " --variants $v"; done`
    java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar \
      -T CombineGVCFs \
      -R ~{refFasta} \
      ${varString} \
      -o ~{outputFileNamePrefix}_jointCalled.g.vcf
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
  }

  output {
    File combinedGVCF = "~{outputFileNamePrefix}_jointCalled.g.vcf"
  }

  parameter_meta {
    vcfs: "Input txt file containing file paths for all VCFs to be combined."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    refFasta: "File path for the human genome reference fasta."
    mem: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
  }

  meta {
    output_meta: {
      combinedGVCF: "Intermediate combined GCVF file"
    }
  }
}

task jointGenotype {
  input {
    String outputFileNamePrefix
    File combinedGCVF = combineGVCF.combinedGVCF
    String? refFasta= "/.mounts/labs/PDE/data/reference/hg19/gatk/2.8/ucsc.hg19.fasta"
    String? modules = "java/8,gatk/3.6-0"
    Int? mem = 32
  }

  command <<<
    java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R ~{refFasta} \
    -V ~{combinedGCVF} \
    -o ~{outputFileNamePrefix}_jointGT.vcf
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
  }

  output {
    File jointGTVCF ="~{outputFileNamePrefix}_jointGT.vcf"
  }

  parameter_meta {
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    combinedGCVF: "Input VCFs generated from combined calling of GCVF."
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
  }

  meta {
    output_meta: {
      jointGTVCF: "Intermedia output joint genotyped VCF."
    }
  }
}

task vqsr {
  input {
    String outputFileNamePrefix
    File jointGTVCF = jointGenotype.jointGTVCF
    String? refFasta= "/.mounts/labs/PDE/data/reference/hg19/gatk/2.8/ucsc.hg19.fasta"
    String?  omniVCF= "/.mounts/labs/PDE/data/gatkAnnotationResources/1000G_omni2.5.hg19.sites.vcf"
    String?  kgSNP= "/.mounts/labs/PDE/data/gatkAnnotationResources/1000g_20100804_chr.vcf"
    String?  hapMapVCF= "/.mounts/labs/PDE/data/gatkAnnotationResources/hapmap_3.3.hg19.sites.vcf"
    String?  dbsnpVCF= "/.mounts/labs/PDE/data/gatkAnnotationResources/dbSNP135_chr.vcf"
    String?  kgIndels= "/.mounts/labs/TGL/gsi/databases/1000G_phase1.indels.b37.vcf_chr.vcf"
    String?  millsIndelsVCF= "/.mounts/labs/TGL/gsi/databases/Mills_and_1000G_gold_standard.indels.b37.vcf_chr.vcf"
    String? tmp="tmp"
    String? modules = "java/8,gatk/3.6-0"
    Int? mem = 32
  }

  command <<<
    # make tmp FOLDER
    mkdir -p ~{tmp}
    # snp model
    java -jar \${GATK_ROOT}/GenomeAnalysisTK.jar \
      -T VariantRecalibrator \
      -R ~{refFasta} \
      -input ~{jointGTVCF} \
      -U ALLOW_SEQ_DICT_INCOMPATIBILITY \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ~{hapMapVCF} \
      --resource:omni,known=false,training=true,truth=true,prior=12.0 ~{omniVCF} \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 ~{kgSNP} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ~{dbsnpVCF} \
      -an DP \
      -an QD \
      -an FS \
      -an SOR \
      -an MQ \
      -an MQRankSum \
      -an ReadPosRankSum \
      -an InbreedingCoeff \
      -mode SNP \
      -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
      -recalFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_SNP.recal \
      -tranchesFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_SNP.tranches \
      -rscriptFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_SNP_plots.R

      # apply snp model
      java -jar \${GATK_ROOT}/GenomeAnalysisTK.jar \
        -T ApplyRecalibration \
        -R ~{refFasta} \
        -input ~{jointGTVCF} \
        -mode SNP \
        --ts_filter_level 99.0 \
        -recalFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_SNP.recal \
        -tranchesFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_SNP.tranches \
        -o ~{outputFileNamePrefix}.recal.snp.indel.vcf

      # build indel model
      java -jar \${GATK_ROOT}/GenomeAnalysisTK.jar \
        -T VariantRecalibrator \
        -R ~{refFasta} \
        -input ~{outputFileNamePrefix}.recal.snp.indel.vcf \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ~{millsIndelsVCF} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ~{kgIndels} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ~{dbsnpVCF} \
        -U ALLOW_SEQ_DICT_INCOMPATIBILITY \
        -an QD \
        -an DP \
        -an FS \
        -an SOR \
        -an MQRankSum \
        -an ReadPosRankSum \
        -an InbreedingCoeff \
        -mode INDEL \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        --maxGaussians 4 \
        -recalFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_INDEL.recal \
        -tranchesFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_INDEL.tranches \
        -rscriptFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_INDEL_plots.R

      # apply indel model
      java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar \
        -T ApplyRecalibration \
        -R ~{refFasta} \
        -input ~{outputFileNamePrefix}.recal.snp.indel.vcf \
        -mode INDEL \
        --ts_filter_level 99.0 \
        -recalFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_INDEL.recal \
        -tranchesFile ~{tmp}/~{outputFileNamePrefix}.recalibrate_INDEL.tranches \
        -o ~{outputFileNamePrefix}.recal.snp.indel.vqsr.vcf
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
  }

  output {
    File vqsrVCF = "~{outputFileNamePrefix}.recal.snp.indel.vqsr.vcf"
    File snpTranche = "~{tmp}/~{outputFileNamePrefix}.recalibrate_SNP.tranches"
    File snpRecal = "~{tmp}/~{outputFileNamePrefix}.recalibrate_ISNP.recal"
    File indelTranche = "~{tmp}/~{outputFileNamePrefix}.recalibrate_INDEL.tranches"
    File indelRecal = "~{tmp}/~{outputFileNamePrefix}.recalibrate_INDEL.recal"
    File snpIndelRecalVcf = "~{outputFileNamePrefix}.recal.snp.indel.vcf"

  }

  parameter_meta {
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    jointGTVCF: "Input VCFs generated from combined calling of GCVF."
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
  }

  meta {
    output_meta: {
      vqsrVCF: "Output VQSR VCF."
      snpTranche: "SNP tranches generated from SNP modelling."
      snpRecal: "SNP recaliberation scores."
      indelTranche: "INDEL tranches generated from INDEL modelling."
      indelRecal: "INDEL recaliberation scores."
      snpIndelRecalVcf: "Intermediate SNP INDEL recaliberation vcf."
    }
  }
}


task computeVariantMetrics {
  input {
    String outputFileNamePrefix
    File vqsrVCF = vqsr.vqsrVCF
    String? dbsnpVCF = "/.mounts/labs/PDE/data/gatkAnnotationResources/dbSNP135_chr.vcf"
    String? modules = "picard/2.19.2"
    Int? mem = 32
  }

  command <<<
    # set up VEP path
    java -jar ${PICARD_ROOT}/picard.jar CollectVariantCallingMetrics \
     INPUT=~{vqsrVCF} \
     OUTPUT=~{outputFileNamePrefix}.recal.snp.indel.vqsr.variant.metrics \
     DBSNP=~{dbsnpVCF}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
  }

  output {
    File variantMetrics = "~{outputFileNamePrefix}.recal.snp.indel.vqsr.variant.metrics"
  }

  parameter_meta {
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    vqsrVCF: "Input VCFs generated from VQSR."
    # modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
  }

  meta {
    output_meta: {
      variantMetrics: "Output VCF metrics."
    }
  }
}

task annotateVCF {
  input {
    String outputFileNamePrefix
    File vqsrVCF = vqsr.vqsrVCF
    String? dbNSFP = "/.mounts/labs/TGL/gsi/databases/vep_plugin_data/dbNSFP_hg19.gz"
    String? refFasta = "/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19/Genomic/references/fasta/hg19.fa"
    String? lofScore = "/.mounts/labs/TGL/gsi/databases/vep_plugin_data/LoFtool_scores.txt"
    # String? modules = "vep/92.0,perl/5.28"
    String? modules = "tabix/1.9"

    Int? mem = 32
  }

  command <<<
    # set up VEP path
    # to be replaced by modulator
    VEP_ROOT=/oicr/local/analysis/sw/vep/vep92
    PERL_ROOT=/oicr/local/analysis/sw/perl/perl-5.22.2-tgl/
    VCF2MAF=/.mounts/labs/PDE/Modules/sw/vcf2maf/mskcc-vcf2maf-decbf60
    # to be replaced by modulator
    export LD_LIBRARY_PATH=${PERL_ROOT}/lib:$LD_LIBRARY_PATH;
    export PERL5LIB=${{PERL_ROOT}/lib:$PERL5LIB
    export PATH=${PERL_ROOT}/bin:$PATH
    # vep/92"
    export PATH=~${VEP_ROOT}:$PATH
    export PATH=${VEP_ROOT}/htslib:$PATH
    export PATH=${VEP_ROOT}/samtools/bin:$PATH
    export PERL5LIB=${VEP_ROOT}:$PATH
    export VEP_PATH=${VEP_ROOT};
    export VEP_DATA=${VEP_ROOT}/.cache;
    # vcf2maf
    export PATH=${VCF2MAF}:$PATH;


    # gzip input VCFs
    cat ~{vqsrVCF} | grep -v END | grep -v GVCFBlock | sed 's/,<NON_REF>//g' | bgzip -c > ~{vqsrVCF}.gz
    tabix -p vcf ~{vqsrVCF}.gz

    # annotate
    ${VEP_ROOT}/vep -i ~{vqsrVCF}.gz \
      --offline --dir_cache $VEP_DATA \
      --assembly GRCh37 \
      --database $VEP_DATA \
      --force_overwrite \
      --fasta ${refFasta} \
      --variant_class \
      --hgvs \
      --symbol \
      --canonical \
      --check_existing \
      --humdiv \
      --sift b \
      --polyphen b \
      --plugin dbNSFP,~{dbNSFP},genename,clinvar_golden_stars,clinvar_clnsig,1000Gp3_AF,ExAC_AF,gnomAD_exomes_AF,gnomAD_genomes_AF,SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred,FATHMM_pred,REVEL_score,CADD_phred,GERP++_RS \
      --plugin LoFtool,~{lofScore} \
      --vcf -o ~{outputFileNamePrefix}.recal.snp.indel.vqsr.annotated.vep.vcf

      cat ~{outputFileNamePrefix}.recal.snp.indel.vqsr.annotated.vep.vcf | bgzip -c > ~{outputFileNamePrefix}.recal.snp.indel.vqsr.annotated.vep.vcf.gz
      tabix -p vcf ~{outputFileNamePrefix}.recal.snp.indel.vqsr.annotated.vep.vcf.gz
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
  }

  output {
    File annotatedVCF = "~{outputFileNamePrefix}.recal.snp.indel.vqsr.annotated.vep.vcf"
  }

  parameter_meta {
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    vqsrVCF: "Input VCFs generated from VQSR."
    # modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    dbNSFP: "Path to dbNSFP file."
    refFasta: "Path to the reference fasta file."
    lofScore: "Path to lofScores txt file."
  }

  meta {
    output_meta: {
      annotatedVCF: "Output annotated VCF."
    }
  }
}
