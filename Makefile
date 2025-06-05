# Makefile for QC process

# Variables
SRR = KWNR_S148_L002_R1_001
ACC = Culex_tarsalis
REF = refs/$(ACC).fa
REF1 = refs/refs_culex_tarsalis/$(ACC).fa
GFF = refs/$(ACC).gff
GTF = refs/$(ACC).gtf
R1 = trim/$(SRR)_trimmed.fastq.gz
R2 = trim2/$(SRR)_trimmed.fastq.gz
BAM = bam/$(SRR).bam
BAM2 = bam2/$(SRR).bam
NA1 = NA_Cu/$(SRR)_unmapped.fq.gz
CONT200 = contigs200.fa


###Important to define this for the sRNA profiling
SAM_FILE = sam200/$(SRR).sam
BAM_FILE = bam200/$(SRR).bam

#######Name of the samples for contig assembly


V1 = KWNR_S148_unmapped_trimmed.fq.gz
V2 = P10_S8_unmapped_trimmed.fq.gz
V3 = P11_S13_unmapped_trimmed.fq.gz
V4 = P12_S9_unmapped_trimmed.fq.gz
V5 = P13_S14_unmapped_trimmed.fq.gz
V6 = P14_S10_unmapped_trimmed.fq.gz
V7 = P15_S15_unmapped_trimmed.fq.gz
V8 = P16_S16_unmapped_trimmed.fq.gz
V9 = P17_S17_unmapped_trimmed.fq.gz
V10 = P1_S45_unmapped_trimmed.fq.gz
V11 = P2_S46_unmapped_trimmed.fq.gz
V12 = P3_S47_unmapped_trimmed.fq.gz
V13 = P4_S4_unmapped_trimmed.fq.gz
V14 = P5_S11_unmapped_trimmed.fq.gz
V15 = P6_S5_unmapped_trimmed.fq.gz
V16 = P7_S6_unmapped_trimmed.fq.gz
V17 = P8_S7_unmapped_trimmed.fq.gz
V18 = P9_S12_unmapped_trimmed.fq.gz



#For the run of the reads that did not align to the reference genome



# Targets and rules
all: folders genome index trim align 
init: folders index 
run: trim align report 


#make all the folders
folders: 
	mkdir -p refs
	mkdir -p trim
	mkdir -p trim2
	mkdir -p reads
	mkdir -p bam
	mkdir -p bam2
	mkdir -p vcf
	mkdir -p merged
	mkdir -p fastqc
	mkdir -p qc_trim
	mkdir -p qc_trim2
	mkdir -p NA_Cu
	mkdir -p NA_Cu2
	mkdir -p NA_ba
	mkdir -p contigs
	mkdir -p contigs_trial
	mkdir -p velvet_output
	mkdir -p sam50
	mkdir -p sam200
	mkdir -p combined_output
	mkdir -p extended_contigs

#quality check of reads with fastqc
fastqc: folders
	fastqc -o fastqc/ reads/$(SRR).fastq.gz

#report of the quality check from all samples
report: fastqc
	multiqc fastqc/ -o fastqc/

###NOTE
#trimming and qc cutadapt. Trims for general use. Removed reads with less than 15 bp. 
trim: folders	
	cutadapt -j 8 --minimum-length=15 -e 0.1 --trimmed-only -a AGATCGGAAGAGC -o trim/$(SRR)_trimmed.fastq.gz reads/$(SRR).fastq.gz


#qc of the trimmed reads. Use R2 only if you have paired-end reads 
qc_trim: folders
	fastqc -o qc_trim/ $(R1)
#	fastqc -o qc_trim2/ $(R2)

#report of the quality check from trimmed reads of all samples
report_trim: qc_trim
	multiqc qc_trim/ -o qc_trim/
#	multiqc qc_trim2/ -o qc_trim2/

# build the index for the reference genome

index:
	@echo "# bowtie2-build index for the reference genome $(REF)"
	bowtie2-build $(REF) $(REF)
	@echo "# Build index for the reference genome $(REF)"
	samtools faidx $(REF)

##Align reads to the reference genome
align:
	bowtie2 -t -p 8 -q --mm -N 0 -x $(REF1) -U $(R1) | samtools view -bS - > $(BAM)
	samtools view -b -f 4 $(BAM) > bam/$(SRR)_unmapped.bam
	samtools fastq bam/$(SRR)_unmapped.bam > NA_Cu/$(SRR)_unmapped.fq
	gzip NA_Cu/$(SRR)_unmapped.fq

#####################Use kraken to map reads to bacteria genome. Better to follow the instructions in the Readme file
bacteria_download:
	kraken2-build --download-taxonomy --db kraken/bact_database
	kraken2-build --download-library bacteria --db kraken/bact_database

Kraken2_build:
	kraken2-build --build --db refs/bact_database --threads 13 --kmer-len 17 --minimizer-len 12 --minimizer-spaces 3

Kraken2:
	kraken2 --db refs/bact_database --threads 32 --unclassified-out NA_ba/$(SRR)_unmapped.fq --report NA_ba/$(SRR)_report.txt --gzip-compressed $(NA1)
	gzip NA_ba/$(SRR)_unmapped.fq

### Contig assembly with velvet and metaspades
#Run velvet, then filter the reads to only those between 20 and 30 bp. Then run velvet with the filtered reads.
velvet:
	zcat NA_ba/$(SRR)_unmapped.fq.gz | awk '{if(NR%4==1) id=$$0; else if(NR%4==2) {seq=$$0; l=length($$0)} else if(NR%4==3) plus=$$0; else if(NR%4==0) {q=$$0; if(l>=20 && l<=30) {print id; print seq; print plus; print q}}}' > NA_ba/$(SRR)_filtered_reads.fq
	mkdir -p velvet_output/$(SRR)
	velveth velvet_output/$(SRR) 15 -fastq -short NA_ba/$(SRR)_filtered_reads.fq
	velvetg velvet_output/$(SRR) -min_contig_lgth 200

###Run metaspades in the NA_ba FOLDER, then merge the contigs of metaspades and velvet. To finally filter out redundant contigs wit CD-HIT
###P9, P11 and P13 only ran with -k 15
spades:
	mkdir -p ../metaspades_output/${SRR}
	spades.py -s $(SRR)_filtered_reads.fq \
	-o ../metaspades_output/${SRR} \
	-t 45 \
	--only-assembler \
	-m 64 \
	--phred-offset 33 \
	-k 15,17 \


###Use cdhit environment (It does not work in the bioinfo environment)
combine_contigs:
	cat velvet_output/${SRR}/contigs.fa metaspades_output/${SRR}/contigs.fasta > combined_output/${SRR}_merged_contigs.fasta
	cd-hit-est -i combined_output/${SRR}_merged_contigs.fasta -o combined_output/${SRR}_contigs.fasta -c 0.95 -n 10

#####Using bioinfo environment, helpts to split the contigs into two files, one with contigs between 50 and 200 bp, and another with contigs larger than 200 bp
split:
	seqtk seq -L 200 combined_output/${SRR}_contigs.fasta > contigs/${SRR}_200.fasta


#######Crate the bowtie index for the contigs, and align the reads to the contigs. Sam files will be necessary for the creation of the sRNA profiles
profile200:
	bowtie2-build renamed_contigs/$(SRR)_renamed.fasta renamed_contigs/$(SRR)_renamed.fasta
	bowtie2 -x renamed_contigs/$(SRR)_renamed.fasta -U NA_ba/$(SRR)_unmapped.fq.gz -S sam200/$(SRR).sam

########split sam files into smaller sam files, one for each contig.

split_sam200:
	mkdir -p sam200/$(SRR)
	samtools view -H sam200/$(SRR).sam header.sam  # Extract header to a separate file
	samtools view -F 4 sam200/$(SRR).sam | awk -v OFS="\t" '{print $$0 > "sam200/$(SRR)/" $$3 ".sam"}'


######get the bam files from original sam files, useful for the coverage plotting
split_bam:
	mkdir -p bam200
	mkdir -p coverage200/$(SRR)
	mkdir -p coverage200/$(SRR)/contigs
	samtools view -bS $(SAM_FILE)> $(BAM_FILE)
	samtools sort $(BAM_FILE) -o $(BAM_FILE)
	samtools index $(BAM_FILE)



########Only use this in the stats folder. It allows us to plot the sRNA profile plots from the contigs. plotGeralDistributionPerBaseByReads.pl comes from  https://github.com/ericgdp/sRNA-virome
stats200:
	mkdir -p sam200
	mkdir -p sam200/$(SRR)
	perl plotGeralDistributionPerBaseByReads.pl -sam ../sam200/$(SRR)/$(CONT200) -s 15 -e 30 -p sam200/$(SRR)/$(CONT200) --plot [--keep]

####################Map samples to file of curated contigs!!!
##########Need to run the following command to create the index for the contigs in the general or root directory
#bowtie2-build non_redundant_contigs.fasta non_redundant_contigs.fasta
##########Coocurrence was used for ultimately producing the contig coocurrence heatmaps
coocurrence:
	mkdir -p co_sam
	mkdir -p co_sam/co_bam
	bowtie2 -N 1 -x non_redundant_contigs.fasta -U NA_ba/$(SRR)_unmapped.fq.gz -S co_sam/$(SRR).sam
	samtools view -bS co_sam/$(SRR).sam > co_sam/co_bam/$(SRR).bam
	samtools sort co_sam/co_bam/$(SRR).bam -o co_sam/co_bam/$(SRR).bam
	samtools index co_sam/co_bam/$(SRR).bam

#######Coocurrence of the reference viruses genomes. This is the one that will be used for the virus coocurrence heatmaps.
coocurrence_final:
	mkdir -p co_sam2
	mkdir -p co_sam2/co_bam
	bowtie2 -N 1 -x refs/reference_viruses.fa -U NA_ba/$(SRR)_unmapped.fq.gz -S co_sam2/$(SRR).sam
	samtools view -bS co_sam2/$(SRR).sam > co_sam2/co_bam/$(SRR).bam
	samtools sort co_sam2/co_bam/$(SRR).bam -o co_sam2/co_bam/$(SRR).bam
	samtools index co_sam2/co_bam/$(SRR).bam

###########################################
###########################################
####Contig extension using the trusted contigs and clusteres found through the clustering of the contigs. Cluster files are in the folder

cluster_all: cluster_1 cluster_2 cluster_3 cluster_4 cluster_5 cluster_6 cluster_7 cluster_8

cluster_1:
	mkdir -p ../extended_contigs/cluster_1
	spades.py \
		-s P1_S45_L002_R1_001_filtered_reads.fq \
		-s P4_S4_L002_R1_001_filtered_reads.fq \
		-s P3_S47_L002_R1_001_filtered_reads.fq \
		-s P2_S46_L002_R1_001_filtered_reads.fq \
		-s P6_S5_L002_R1_001_filtered_reads.fq \
		-s P13_S14_L002_R1_001_filtered_reads.fq \
		-s P15_S15_L002_R1_001_filtered_reads.fq \
		-s P17_S17_L002_R1_001_filtered_reads.fq \
		-s P7_S6_L002_R1_001_filtered_reads.fq \
		-o ../extended_contigs/cluster_1 \
		-t 4 \
		--only-assembler \
		-m 32 \
		--phred-offset 33 \
		-k 15,17 \
		--careful \
		--trusted-contigs cluster_1.fasta

cluster_2:
	mkdir -p ../extended_contigs/cluster_2
	spades.py \
		-s P1_S45_L002_R1_001_filtered_reads.fq \
		-s P4_S4_L002_R1_001_filtered_reads.fq \
		-s P13_S14_L002_R1_001_filtered_reads.fq \
		-o ../extended_contigs/cluster_2 \
		-t 4 \
		--only-assembler \
		-m 32 \
		--phred-offset 33 \
		-k 15,17 \
		--careful \
		--trusted-contigs cluster_2.fasta

cluster_3:
	mkdir -p ../extended_contigs/cluster_3
	spades.py \
		-s P1_S45_L002_R1_001_filtered_reads.fq \
		-s P4_S4_L002_R1_001_filtered_reads.fq \
		-s P3_S47_L002_R1_001_filtered_reads.fq \
		-s P2_S46_L002_R1_001_filtered_reads.fq \
		-s P6_S5_L002_R1_001_filtered_reads.fq \
		-s P8_S7_L002_R1_001_filtered_reads.fq \
		-s P5_S11_L002_R1_001_filtered_reads.fq \
		-s P16_S16_L002_R1_001_filtered_reads.fq \
		-s P7_S6_L002_R1_001_filtered_reads.fq \
		-s P14_S10_L002_R1_001_filtered_reads.fq \
		-o ../extended_contigs/cluster_3 \
		-t 4 \
		--only-assembler \
		-m 32 \
		--phred-offset 33 \
		-k 15,17 \
		--careful \
		--trusted-contigs cluster_3.fasta

cluster_4:
	mkdir -p ../extended_contigs/cluster_4
	spades.py \
		-s P1_S45_L002_R1_001_filtered_reads.fq \
		-s P2_S46_L002_R1_001_filtered_reads.fq \
		-s P3_S47_L002_R1_001_filtered_reads.fq \
		-s P4_S4_L002_R1_001_filtered_reads.fq \
		-s P5_S11_L002_R1_001_filtered_reads.fq \
		-s P6_S5_L002_R1_001_filtered_reads.fq \
		-s P7_S6_L002_R1_001_filtered_reads.fq \
		-s P8_S7_L002_R1_001_filtered_reads.fq \
		-s P9_S12_L002_R1_001_filtered_reads.fq \
		-s P11_S13_L002_R1_001_filtered_reads.fq \
		-s P12_S9_L002_R1_001_filtered_reads.fq \
		-s P13_S14_L002_R1_001_filtered_reads.fq \
		-s P14_S10_L002_R1_001_filtered_reads.fq \
		-s P15_S15_L002_R1_001_filtered_reads.fq \
		-s P16_S16_L002_R1_001_filtered_reads.fq \
		-s P17_S17_L002_R1_001_filtered_reads.fq \
		-o ../extended_contigs/cluster_4 \
		-t 4 \
		--only-assembler \
		-m 32 \
		--phred-offset 33 \
		-k 15,17 \
		--careful \
		--trusted-contigs cluster_4.fasta

cluster_5:
	mkdir -p ../extended_contigs/cluster_5
	spades.py \
		-s P1_S45_L002_R1_001_filtered_reads.fq \
		-s P4_S4_L002_R1_001_filtered_reads.fq \
		-s P11_S13_L002_R1_001_filtered_reads.fq \
		-s P3_S47_L002_R1_001_filtered_reads.fq \
		-s P2_S46_L002_R1_001_filtered_reads.fq \
		-s P6_S5_L002_R1_001_filtered_reads.fq \
		-s P13_S14_L002_R1_001_filtered_reads.fq \
		-s P9_S12_L002_R1_001_filtered_reads.fq \
		-s P8_S7_L002_R1_001_filtered_reads.fq \
		-o ../extended_contigs/cluster_5 \
		-t 4 \
		--only-assembler \
		-m 32 \
		--phred-offset 33 \
		-k 15,17 \
		--careful \
		--trusted-contigs cluster_5.fasta

cluster_6:
	mkdir -p ../extended_contigs/cluster_6
	spades.py \
		-s P1_S45_L002_R1_001_filtered_reads.fq \
		-s P4_S4_L002_R1_001_filtered_reads.fq \
		-s P11_S13_L002_R1_001_filtered_reads.fq \
		-s P3_S47_L002_R1_001_filtered_reads.fq \
		-s P2_S46_L002_R1_001_filtered_reads.fq \
		-s P6_S5_L002_R1_001_filtered_reads.fq \
		-s P17_S17_L002_R1_001_filtered_reads.fq \
		-s P7_S6_L002_R1_001_filtered_reads.fq \
		-s P14_S10_L002_R1_001_filtered_reads.fq \
		-s P12_S9_L002_R1_001_filtered_reads.fq \
		-o ../extended_contigs/cluster_6 \
		-t 4 \
		--only-assembler \
		-m 32 \
		--phred-offset 33 \
		-k 15,17 \
		--careful \
		--trusted-contigs cluster_6.fasta

cluster_7:
	mkdir -p ../extended_contigs/cluster_7
	spades.py \
		-s P1_S45_L002_R1_001_filtered_reads.fq \
		-s P4_S4_L002_R1_001_filtered_reads.fq \
		-s P11_S13_L002_R1_001_filtered_reads.fq \
		-s P3_S47_L002_R1_001_filtered_reads.fq \
		-s P2_S46_L002_R1_001_filtered_reads.fq \
		-s P6_S5_L002_R1_001_filtered_reads.fq \
		-s P7_S6_L002_R1_001_filtered_reads.fq \
		-s P14_S10_L002_R1_001_filtered_reads.fq \
		-o ../extended_contigs/cluster_7 \
		-t 4 \
		--only-assembler \
		-m 32 \
		--phred-offset 33 \
		-k 15,17 \
		--careful \
		--trusted-contigs cluster_7.fasta

cluster_8:
	mkdir -p ../extended_contigs/cluster_8
	spades.py \
		-s P1_S45_L002_R1_001_filtered_reads.fq \
		-s P4_S4_L002_R1_001_filtered_reads.fq \
		-s P11_S13_L002_R1_001_filtered_reads.fq \
		-s P3_S47_L002_R1_001_filtered_reads.fq \
		-s P2_S46_L002_R1_001_filtered_reads.fq \
		-s P6_S5_L002_R1_001_filtered_reads.fq \
		-s P13_S14_L002_R1_001_filtered_reads.fq \
		-s P9_S12_L002_R1_001_filtered_reads.fq \
		-s P5_S11_L002_R1_001_filtered_reads.fq \
		-s P7_S6_L002_R1_001_filtered_reads.fq \
		-s P14_S10_L002_R1_001_filtered_reads.fq \
		-s P12_S9_L002_R1_001_filtered_reads.fq \
		-o ../extended_contigs/cluster_8 \
		-t 4 \
		--only-assembler \
		-m 32 \
		--phred-offset 33 \
		-k 15,17 \
		--careful \
		--trusted-contigs cluster_8.fasta


#####################ONLY used for sRNA profiling. Just align everything to a file containing all the reference genomes. Allow 1 mismatch (-N 1)
#and provide all possible alignments

coocurrence_all:
	mkdir -p viruses
	bowtie2 -t -p 8 -q --mm -N 1 -k 10 -x $(REF) -U NA_ba/P10_S8_L002_R1_001_unmapped.fq.gz,NA_ba/P11_S13_L002_R1_001_unmapped.fq.gz,NA_ba/P12_S9_L002_R1_001_unmapped.fq.gz,NA_ba/P13_S14_L002_R1_001_unmapped.fq.gz,NA_ba/P14_S10_L002_R1_001_unmapped.fq.gz,NA_ba/P15_S15_L002_R1_001_unmapped.fq.gz,NA_ba/P16_S16_L002_R1_001_unmapped.fq.gz,NA_ba/P17_S17_L002_R1_001_unmapped.fq.gz,NA_ba/P1_S45_L002_R1_001_unmapped.fq.gz,NA_ba/P2_S46_L002_R1_001_unmapped.fq.gz,NA_ba/P3_S47_L002_R1_001_unmapped.fq.gz,NA_ba/P4_S4_L002_R1_001_unmapped.fq.gz,NA_ba/P5_S11_L002_R1_001_unmapped.fq.gz,NA_ba/P6_S5_L002_R1_001_unmapped.fq.gz,NA_ba/P7_S6_L002_R1_001_unmapped.fq.gz,NA_ba/P8_S7_L002_R1_001_unmapped.fq.gz,NA_ba/P9_S12_L002_R1_001_unmapped.fq.gz -S viruses/all.sam
	samtools view -bS viruses/all.sam > viruses/all.bam
	samtools sort viruses/all.bam -o viruses/all.bam
	samtools index viruses/all.bam

########Only use this in the stats folder. This will create the sRNA profile plots from the viruses.
stats_virus:
	mkdir -p ../viruses/$(SRR)
	perl plotGeralDistributionPerBaseByReads.pl -sam ../viruses/$(SRR).sam -s 15 -e 30 -p ../viruses/$(SRR)/$(SRR) --plot [--keep]

#########Same code, but using it in the overal alignment files (all the libraries together against the reference genomes, then splited into
#individual sam files per viral genome

####Use this
stats_virus2:
	mkdir -p ../viruses/split/$(SRR)
	perl plotGeralDistributionPerBaseByReads.pl -sam ../viruses/split/$(SRR).sam -s 15 -e 30 -p ../viruses/split/$(SRR)/$(SRR) --plot [--keep]

viruses_bam:
	mkdir -p viruses/bam
	samtools view -bS viruses/$(SRR).sam > viruses/bam/$(SRR).bam
	samtools sort viruses/bam/$(SRR).bam -o viruses/bam/$(SRR).bam
	samtools index viruses/bam/$(SRR).bam

###########################Same than before, but in the folder viruses/split
viruses_bam2:
	mkdir -p viruses/split/bam
	samtools view -bS viruses/split/$(SRR).sam > viruses/split/bam/$(SRR).bam
	samtools sort viruses/split/bam/$(SRR).bam -o viruses/split/bam/$(SRR).bam
	samtools index viruses/split/bam/$(SRR).bam

