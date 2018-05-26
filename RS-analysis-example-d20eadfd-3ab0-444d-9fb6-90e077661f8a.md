# RNAseq analysis example

An outline of my code used to analyze Timelapse-seq data

    # download data from ftp server
    mkdir timelapse-seq
    cd timelapse-seq
    wget -m --user=gslftp --password=gsl23ftp!ftp://gslserver.qb3.berkeley.edu/180427_50SR_HS4KA/*
    cd gslserver.qb3.berkeley.edu/L2345678/
    mv Ingolia-RSLiB/ ~/timelapse-seq/
    cd ~/timelapse-seq/
    rm -r gslserver.qb3.berkeley.edu/
    mv Ingolia-RSLiB/ RawFastq
    
    # removing adaptor sequences from library, then splitting by barcode
    cd ~/timelapse-seq
    ./adapter_removal.sh
    
    # 
    cd ~/timelapse-seq/Alignment
    ln -s ~/timelapse-seq/RawFastq/*_clipped.fastq .
    

## Scripts used to trim adapters and map reads

- **Adapter_removal.sh**

  I made a script to run more easily going forward called `adapter_removal.sh`. I'm placing the Ribosome Profiling libraries and the RNAseq libraries into bash arrays to loop over using different parameters at the `fastx_clipper` (-c) stage. 

      #!/bin/bash
      
      # prep dir
      mkdir -p RawFastq
      mkdir -p Alignment
      mkdir -p Analysis
      
      
      # enter directory name of Ribosome Profilng & RNAseq samples
      cd ./RawFastq
      RPSEQ=(RP*.gz)
      RNASEQ=(RS*.gz)
      DIR='/mnt/ingolialab/apadron/timelapse-seq/RawFastq/'
      ILLUMINA='AGATCGGAAGAGCACACGTCTGAA'
      
      for i in "${RNASEQ[@]}"; do
      			nohup zcat $i | fastx_clipper -Q33 -a ${ILLUMINA} \
      							-v -o ${i}_clipped.fastq &
      done
      
      for i in "${RPSEQ[@]}"; do
      	nohup zcat $i | fastq_illumina_filter --keep N -v | \
      	fastx_clipper -Q33 -a ${ILLUMINA} \
      	-c -v -o ${i}_clipped.fastq &
      done 
      
      # stats from removing barcode
      for i in *_clipped.fastq; do
      		zcat $i | wc -l
      done
      
      wc -l ${i}_clipped.fastq

- **RP_barcode_split.sh**

  This requires having a **pool.csv** file with the names of each file and the corresponding barcode sequences. 

  [pool.csv](https://www.notion.so/file/https%3A%2F%2Fs3-us-west-2.amazonaws.com%2Fsecure.notion-static.com%2F7f147045-bd68-46fe-8d14-7db4fa2d351e%2Fpool.csv)

      #!/bin/bash
      
      cd Alignment
      RPDIRNAME='RP'
      RPDIR='/mnt/ingolialab/apadron/timelapse-seq/Alignment/'
      
      RPLIB=(RP*_clipped.fastq)
      rRNADIR='/mnt/ingolialab/apadron/Genomes/Homo_sapiens/rRNARef/'
      
      # splitting reads along with the barcode
      for i in "${RPLIB[@]}"; do
      nohup fastx-split -o ${RPDIRNAME} -p NN -x NNNNNIIIII --min-insert=10 -s pool.csv ${i} &
      done
      
      # bowties index
      #cd /mnt/ingolialab/apadron/Genomes/Homo_sapiens/rRNARef/
      #bowtie2-build -f RefSeqrRNA.fa RefSeqrRNA
      
      # stats on RP barcode splitting
      cd ${RPDIR}${RPDIRNAME}
      cat fates.txt
      
      
      
      # align reads to human rRNA database and capture un-aligned reads
      cd ${RPDIR}${RPDIRNAME}
      for i in *.fastq; do
      bowtie2 -p 2 --very-sensitive --quiet --un ${i}.norrna.fq -x ${rRNADIR}RefSeqrRNA -U ${i} | rrna-stats -o ${i}.stats --tam --maxread 100 --lenrange 5,100 - &
      done

- **testing_alignment_param.sh**

  **Note: This is the code I proceeded with**

  This script is used to test the `--ignore-quals` option in Hisat2, since the frequency of G>A mutations was low in some of the test samples

      #!/bin/bash
      
      MAINDIR='/mnt/ingolialab/apadron/timelapse-seq'
      cd ${MAINDIR}
      INDEX='/mnt/ingolialab/apadron/Genomes/Homo_sapiens/GRCh38/hisat/GRCh38'
      # define RNAseq dir and put files into array
      RNASEQDIR='/mnt/ingolialab/apadron/timelapse-seq/Alignment/testing_other_alignment_params/'
      RNASEQ=(`find ${RNASEQDIR} -name RS*_clipped.fastq`)
      
      # define RPseq dir and put files into array
      RPSEQDIR='/mnt/ingolialab/apadron/timelapse-seq/Alignment/testing_other_alignment_params'
      RPSEQ=(`find ${RPSEQDIR} -name *.norrna.fq`)
      
      
      # map reads using Hisat2
      # Convert SAM to BAM format
      for i in "${RNASEQ[@]}"; do 
      hisat2 -p 32 --ignore-quals --score-min L,0.0,-0.5 -q -x ${INDEX} -U ${i} -S ${i}.sam
      # Convert SAM to BAM format
      samtools view -bS ${i}.sam > ${i}.bam
      done
      
      
      # map reads using Hisat2
      # Convert SAM to BAM format
      for i in "${RPSEQ[@]}"; do 
      hisat2 -p 32 --ignore-quals --score-min L,0.0,-0.5 -q -x ${INDEX} -U ${i} -S ${i}.sam
      samtools view -bS ${i}.sam > ${i}.bam
      done
      
      # Sort BAM files
      BAMFILES=(`find ./Alignment/testing_other_alignment_params/* -name *.bam`)
      for i in "${BAMFILES[@]}"; do
      	samtools sort $i -o ${i}.sorted.bam &
      done

- **Run FeatureCounts**

  Count the number of reads mapping to each transcript

      #!/bin/bash
      
      # run FeatureCounts on ALL BAM files
      cd ~/timelapse-seq/Analysis
      
      FEATURECOUNTSDIR='/mnt/ingolialab/apadron/subread-1.5.3-Linux-x86_64/bin/'
      GTFDIR='/mnt/ingolialab/apadron/Genomes/Homo_sapiens/GRCh38/'
      
      for i in *.sorted.bam ; do
      ${FEATURECOUNTSDIR}featureCounts -t exon -M --fraction -g gene_id \
      -a ${GTFDIR}gencode.gtf -o ${i}.featureCounts ${i} &
      done
      
      for i in *.sorted-yestl.bam; do
      ${FEATURECOUNTSDIR}featureCounts -t exon -M --fraction -g gene_id \
      -a ${GTFDIR}gencode.gtf -o ${i}.featureCounts ${i} &
      done