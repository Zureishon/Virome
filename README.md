# This serves as the protocol for the methods used in the paper https://doi.org/10.1101/2025.05.23.655811

First, we need to get the reference genome for Culex tarsalis (Referenced in paper), and transform the gff3 file into gff and gtf.
(not specified). and save it into a /refs directory. in my case, the root directory was /Virome

All the necessary codes are either in the Makefile, or in R scripts. Makefile has to be in the root directory, or in a different one if specified.
Here I use pre-existing environments found in the book Biostar https://www.biostarhandbook.com/
Or some custom ones when the apps were not compatible with existing environments, as it was the case for cdhit and velvet

#### Note:
If I have any problem with sbatch scripts, run the dos2unix command

Example:

dos2unix assembly_make.sh


## Qc process of the raw reads (trimmed by novogene)

cat samples.txt | parallel --lb make fastqc SRR={}

### multiqc

make report

## initial trimming

cat samples.txt | parallel --lb make trim SRR={}

## QC of trimmed reads

cat samples.txt | parallel --lb make qc_trim SRR={}

## Run for creating the C. tarsalis genome bowtie index

make index ACC=Culex_tarsalis

## Run for alignment to C. tarsalis reference genome. Unmapped reads will be saved 
make align SRR=KWNR_S148_L002_R1_001  
make align SRR=P10_S8_L002_R1_001  
make align SRR=P11_S13_L002_R1_001  
make align SRR=P12_S9_L002_R1_001  
make align SRR=P13_S14_L002_R1_001   
make align SRR=P14_S10_L002_R1_001  
make align SRR=P15_S15_L002_R1_001  
make align SRR=P16_S16_L002_R1_001  
make align SRR=P17_S17_L002_R1_001  
make align SRR=P1_S45_L002_R1_001  
make align SRR=P2_S46_L002_R1_001  
make align SRR=P3_S47_L002_R1_001  
make align SRR=P4_S4_L002_R1_001  
make align SRR=P5_S11_L002_R1_001  
make align SRR=P6_S5_L002_R1_001  
make align SRR=P7_S6_L002_R1_001  
make align SRR=P8_S7_L002_R1_001  
make align SRR=P9_S12_L002_R1_001  

or  

cat samples.txt | parallel --lb make align SRR={}  

## Create Kraken2 reference database  
make bacteria_download  
make Kraken2_build  

## Then align the remaining reads to bacteria genome using Kraken2  
cat samples.txt | parallel --lb make Kraken2 SRR={}  

## Contig assembly

cat samples.txt | parallel --lb make velvet SRR={}  
cat samples.txt | parallel --lb make spades SRR={}  

##### For samples SRR=P9_S12_L002_R1_001, SRR=P13_S14_L002_R1_001, SRR=P11_S13_L002_R1_001, the contig assembly can only be done with -k 15  

## Now run cdhit and split (To get the 200bp contigs in the contigs folder). It basically pulls the contigs from velvet and SPAdes together.  
cat samples.txt | parallel --lb make combine_contigs SRR={}  
cat samples.txt | parallel --lb make split SRR={}  

### After these steps, I manually blasted each file contigs/${SRR}_200.fasta, downloaded the csv file from NCBI blast, and ran the Rscript "blast_format.R" (located in the contigs folder) for each sample  
### Then put the resulting contigs files (renamed), into the directory renamed_contigs   

## Alignment of reads to the contigs 
cat samples200.txt | parallel --lb make profile200 SRR={}  
cat samples200.txt | parallel --lb make split_sam200 SRR={}  

## Create sRNA profiles for each contig of each sam file. Use in stats environment and directory. The directory must include the R and pearl scripts from Aguiar et. al  
ls ../sam200/KWNR_S148_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=KWNR_S148_L002_R1_001 CONT200={}  
ls ../sam200/P1_S45_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P1_S45_L002_R1_001 CONT200={}  
ls ../sam200/P2_S46_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P2_S46_L002_R1_001 CONT200={}  
ls ../sam200/P3_S47_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P3_S47_L002_R1_001 CONT200={}  
ls ../sam200/P4_S4_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P4_S4_L002_R1_001 CONT200={}  
ls ../sam200/P5_S11_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P5_S11_L002_R1_001 CONT200={}  
ls ../sam200/P6_S5_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P6_S5_L002_R1_001 CONT200={}  
ls ../sam200/P7_S6_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P7_S6_L002_R1_001 CONT200={}  
ls ../sam200/P14_S10_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P14_S10_L002_R1_001 CONT200={}  
ls ../sam200/P15_S15_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P15_S15_L002_R1_001 CONT200={}  
ls ../sam200/P17_S17_L002_R1_001/ | parallel --lb make -f ../Makefile stats200 SRR=P17_S17_L002_R1_001 CONT200={}  

## Here comes the manual curation. I suggest moving all the profiles into a single folder, and double check the virome profiles looking for 21bp peaks.  
In this case, those are the contigs per sample that seemed to have viral origin  

##KWNR =   
None  
##P1 =  
P1_contig006_978_MAG
P1_contig012_422_MAG
P1_contig013_306_MAG
P1_contig014_283_Culex
P1_contig016_276_MAG
P1_contig017_251_MAG
P1_contig019_244_MAG
P1_contig020_242_Culex  
##P2 =  
P2_contig001_434_Partitivirus-like
P2_contig002_477_NA_NA
P2_contig003_328_MAG
P2_contig004_205_MAG
P2_contig005_239_Culex
P2_contig006_1013_NA_NA
P2_contig007_537_NA_NA
P2_contig010_416_Partitivirus-like
P2_contig011_387_NA_NA
P2_contig012_379_MAG
P2_contig014_354_Partitivirus-like
P2_contig015_335_Partitivirus-like
P2_contig017_284_NA_NA
P2_contig019_250_NA_NA
P2_contig021_239_NA_NA  
##P3 =   
P3_contig008_268_NA_NA  
##P4 =  
P4_contig003_268_NA_NA  
##P6 =   
P6_contig003_750_NA_NA  
##P7 =   
P7_contig001_419_NA_NA  
##P17 =   
P17_contig001_449_Culex  

## It is important to create a file called curated_contigs.fasta that includes such contigs   

## Then use this code to decrease redundancy of contigs (use in root directory)  
cd-hit -i curated_contigs.fasta -o non_redundant_contigs.fasta -c 0.90 -aS 0.90  

## Map the remaining reads (mosquito and bacteria free reads) to the contigs for coocurrence analysis   
cat samples.txt | parallel --lb make coocurrence SRR={}  

## Calculate the readcounts  
featureCounts -F GTF -t contig -a non_redundant_contigs.gtf -o co_sam/co_bam/counts.txt co_sam/co_bam/KWNR_S148_L002_R1_001.bam co_sam/co_bam/P14_S10_L002_R1_001.bam co_sam/co_bam/P15_S15_L002_R1_001.bam co_sam/co_bam/P17_S17_L002_R1_001.bam co_sam/co_bam/P1_S45_L002_R1_001.bam co_sam/co_bam/P2_S46_L002_R1_001.bam co_sam/co_bam/P3_S47_L002_R1_001.bam co_sam/co_bam/P4_S4_L002_R1_001.bam co_sam/co_bam/P5_S11_L002_R1_001.bam co_sam/co_bam/P6_S5_L002_R1_001.bam co_sam/co_bam/P7_S6_L002_R1_001.bam co_sam/co_bam/P8_S7_L002_R1_001.bam co_sam/co_bam/P10_S8_L002_R1_001.bam co_sam/co_bam/P11_S13_L002_R1_001.bam co_sam/co_bam/P12_S9_L002_R1_001.bam co_sam/co_bam/P13_S14_L002_R1_001.bam co_sam/co_bam/P16_S16_L002_R1_001.bam co_sam/co_bam/P9_S12_L002_R1_001.bam

### Use the resulting counts.txt file and RPKM.R script to make the heatmaps, and retrieve the contig and reads clusters for contig extension with spades. Then take all the cluster_#.fasta files into the root directory 

## Run the scripts in the makefile for each cluster (modify as needed)  
make cluster_1  
make cluster_2  
make cluster_3  
make cluster_4  
make cluster_5  
make cluster_6  
make cluster_7  
make cluster_8  

## In the root folder use for combining all the extended contigs, and keeping only the ones longer than 200bp, and sort them by size  
cat extended_contigs/*/contigs.fasta > final_contigs/contigs.fasta  
cd-hit-est -i final_contigs/contigs.fasta -o final_contigs/contigs2.fasta -c 0.95 -n 10  
seqtk seq -L 200 final_contigs/contigs2.fasta > final_contigs/contigs_final.fasta  
seqkit sort -l -r final_contigs/contigs_final.fasta > final_contigs/contigs_final_sorted.fasta  

## From here, we repeat the steps for blasting and use of the blast_format.R, but only for the contigs_final_sorted.fasta file  
This step allowed us to match the extended contigs with 7 insect specific viruses. Reference genomes were selected and downloaded according to the blast results.  

MH188052.1,Culex_Bunyavirus_2  
NC_040716.1,Culex_Iflavi-like_virus_4                  same as MH188011.1 (Genbank)  
NC_032231.1,Hubei_mosquito_virus_4                     same as KX883008.1 (Genbank)  
MH188050.1,Partitivirus-like_Culex_mosquito_virus  
MF176248.1,Wuhan_Mosquito_Virus_6  
MW434901.1,Marma_virus  
MK628543.1,Culex_narnavirus_1  

## Add the reference genomes into a fasta file called reference_viruses.fa, and index it   
bowtie2-build refs/reference_viruses.fa refs/reference_viruses.fa  

## Alignment required for sRNA profiling of viruses  
make coocurrence_all ACC=reference_viruses  

## Then remove all secondary and chimeric alignments  
samtools view -h -F 256 -F 2048 viruses/all.sam > viruses/all_filtered.sam  
samtools sort viruses/all_filtered.bam -o viruses/all_filtered.bam  
samtools index viruses/all_filtered.bam  

## Split sam files  
for contig in $(samtools view -H all_filtered.bam | grep '^@SQ' | awk '{print $2}' | cut -d ':' -f 2); do samtools view -h all_filtered.bam $contig > split/${contig}.sam; done
  
## Create the sRNA profiles for each virus. Go to stats folder   

make -f ../Makefile stats_virus2 SRR=Culex_Bunyavirus_2  
make -f ../Makefile stats_virus2 SRR=Culex_Iflavi-like_virus_4  
make -f ../Makefile stats_virus2 SRR=Hubei_mosquito_virus_4  
make -f ../Makefile stats_virus2 SRR=Partitivirus_like_virus_4  
make -f ../Makefile stats_virus2 SRR=Wuhan_Mosquito_Virus_6  
make -f ../Makefile stats_virus2 SRR=Marma_virus  
make -f ../Makefile stats_virus2 SRR=Culex_narnavirus_1  

## Next is the alignment for coocurrence of virus

For this is necessary a refs/reference_viruses.fa file, which is the one containing the reference genomes of the viruses. 

### First create the index for the reference viruses  
bowtie2-build refs/reference_viruses.fa refs/reference_viruses.fa  

### Then align the reads to the reference viruses  
cat samples.txt | parallel --lb make coocurrence_final SRR={}  

### For creating the coverage plots, take the bam and .bai files from the co_sam2/co_bam/ directory, and run the R script "coverage.R"  


