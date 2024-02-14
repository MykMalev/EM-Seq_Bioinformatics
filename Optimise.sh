#Code used in optimisation of ASM pipeline.
#The selected chromosome for testing was chr20. This is because it has a well-defined imprinted region GNAS which serves as a control for Allele-Specific methylation. In addition, it is one of the shortest chromosomes, thus making the preliminary analysis and optimization faster.

#WGS

fastqc WGS.R1
fastqc WGS.R2
trim_galore -o ~/WGS_Trimmed --gzip --paired WGS.R1 WGS.R2
bwa mem -o ~/WGS.Mapped ~HG38.fa WGS.R1 WGS.R2
samtools view -b ~/WGS.mapped > ~/WGS.mapped.bam
samtools sort ~/WGS.sorted.bam
picard MarkDuplicates I=WGS.sorted.bam O=WGS.deduplicated M=WGS.Deduplication_Metrics REMOVE_DUPLICATES=TRUE
gatk HaplotypeCaller -mbq 20 -I ~/WGS.deduplicated -O ~/WGS.SNPs.vcf -R ~/HG38.fa
bcftools view ~/WGS.SNPs.vcf.gz --regions chr20 > WGS.SNPs.chr20.vcf

#Calling SNPs using BS-SNPer
perl ~/BS-Snper.pl ~/chr20.bam --fa ~/HG38_with_Ctrl --output ~/chr20.BS_SNPer.snp.vcf > ~/chr20.out 2> ~/chr20ERR.log &

#Calling SNPs using bis-snp
perl ~/bissnp_easy_usage.pl --nt 8 ~/BisSNP-0.77.jar ~/chr20.bam ~/HG38_with_Ctrl.fa ~/dbsnp/00-common_all.vcf.gz

#### CGmapTools

#ASM with BS-SNPer SNPs
#ASR

cgmaptools asm -o ~/chr20.BS_SNPer.ass -r ~/HG38_with_Ctrl.fa -b ~/chr20.bam -l ~/chr20.BS_SNPer.snp.vcf

#ASS
cgmaptools asm -o ~/chr20.BS_SNPer.ass -m ass -r ~/HG38_with_Ctrl.fa -b ~/chr20.bam -l ~/chr20.BS_SNPer.snp.vcf


#ASM with Bis-SNP SNPs
#ASR
cgmaptools asm -o ~/chr20.bis-snp.asr  -r ~/HG38_with_Ctrl.fa -b ~/chr20.bam -l ~/chr20.bis_snp.snp.vcf

#ASS
cgmaptools asm -o ~/chr20.bis-snp.ass -m ass -r ~/HG38_with_Ctrl.fa -b ~/chr20.bam -l ~/chr20.bis_snp.snp.vcf

#CGmapTools analysis using WGS SNPs called by GATK


#ASR
cgmaptools asm -o ~/CGmapTools/WGS.chr20.asr -r ~/HG38 -b ~/chr20.bam -l ~/WGS.SNPs.chr20.vcf


#ASS

cgmaptools asm -o ~/CGmapTools/WGS.chr20.asr -m ass  -r ~/HG38 -b ~/chr20.bam -l ~/WGS.SNPs.chr20.vcf



#ASS to Bedgraph> 
awk ' {print $1, $6, $9, $10}' OFS='\t' ${ass} > ${ass}.tmp
LC_NUMERIC=C awk ' {if ($3>$4) print $1, $2, $3, $4, $3-$4; else print $1, $2, $3, $4, $4-$3}' OFS='\t' ${ass}.tmp > ${ass}.tmp2
awk ' {print $1, $2, $2+1, $5}' OFS='\t' ${ass}.tmp2 > ${ass}.bedgraph
 
#ASR to bedgraph
awk ' {print $1, $6, $7, $10, $11}' OFS='\t' ${asr} > ${asr}.tmp
LC_NUMERIC=C awk ' {if ($5>$4) print $1, $2, $3, $4,$5, $5-$4; else print $1, $2, $3, $4,$5, $4-$5}' OFS='\t' ${asr}.tmp > ${asr}.tmp2
awk ' {print $1, $2, $3, $6}' OFS='\t' ${asr}.tmp2 > [asr].bedgraph


# DNMTools
~/dnmtools states -z  -o ./chr20.epiread -t 4 -c ~/HG38_with_Ctrl.fa ~/chr20.bam

~/dnmtools allelic -v -o ~/chr20.allelic -c ~/HG38_with_Ctrl.fa ~/chr20.epiread

~/dnmtools amrfinder -o ~/chr20.AMR -c ~/HG38_with_Ctrl.fa -m 10 chr20.epiread



#DAMEFinder (RStudio)


library(DAMEfinder)
library(GenomicRanges)
library(SummarizedExperiment)
DATA_PATH_DIR <- system.file('extdata', package = 'DAMEfinder')
get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)
ref <- get_data_path('Combined_Genome.fa')




#Using WGS SNPs

bamFiles <- get_data_path("chr20.bam")
WGSvcfFile <- get_data_path("WGS.SNPs.chr20.vcf")
sampleNames <- "Chr20"
WGS_chr20_extractbam <- extract_bams(bamFiles, WGSvcfFile, sampleNames, coverage = 6, ref, cores = 32)

#Extract ASM scores
WGS_chr20_derASM <- calc_derivedasm(WGS_chr20_extractbam, cores = 32, verbose = TRUE)
WGS_chr20_asm <- assay(WGS_chr20_derASM, "der.ASM" )
write.table(WGS_chr20_asm, file = "WGS_chr20asm.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#Extract SNPs linked
WGS_chr20_snp <- assay(WGS_chr20_derASM, "snp.table")
write.table(WGS_chr20_snp, file = "WGS_chr20_snp.txt" , sep = "\t",
            row.names = TRUE, col.names = NA)
#Extract zASM score
WGS_chr20_zASM <- assay(WGS_chr20_derASM, "zASM")
write.table(WGS_chr20_zASM, file = "WGS_chr20_zASM.txt" , sep = "\t",
            row.names = TRUE, col.names = NA)


#Using BS-SNPer SNPs

BSSNPervcfFile <- get_data_path("chr20.BS_SNPer.snp.vcf")
sampleNames <- "Chr20"
BSSNPer_chr20_extractbam <- extract_bams(bamFiles, BSSNPervcfFile, sampleNames, coverage = 6, ref, cores = 32)

#Extract ASM scores
BSSNPer_chr20_derASM <- calc_derivedasm(BSSNPer_chr20_extractbam, cores = 32, verbose = TRUE)
BSSNPer_chr20_asm <- assay(BSSNPer_chr20_derASM, "der.ASM" )
write.table(BSSNPer_chr20_asm, file = "BSSNPer_chr20asm.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#Extract SNPs linked
BSSNPer_chr20_snp <- assay(BSSNPer_chr20_derASM, "snp.table")
write.table(BSSNPer_chr20_snp, file = "BSSNPer_chr20_snp.txt" , sep = "\t",
            row.names = TRUE, col.names = NA)
#Extract zASM score
BSSNPer_chr20_zASM <- assay(BSSNPer_chr20_derASM, "zASM")
write.table(BSSNPer_chr20_zASM, file = "BSSNPer_chr20_zASM.txt" , sep = "\t",
            row.names = TRUE, col.names = NA)



#Using Bis-snp SNPs
BisSNPvcfFile <- get_data_path("chr20.bis_snp.snp.vcf")
sampleNames <- "Chr20"
BisSNP_chr20_extractbam <- extract_bams(bamFiles, BisSNPvcfFile, sampleNames, coverage = 6, ref, cores = 32)

#Extract ASM scores
BisSNP_chr20_derASM <- calc_derivedasm(BisSNP_chr20_extractbam, cores = 32, verbose = TRUE)
BisSNP_chr20_asm <- assay(BisSNP_chr20_derASM, "der.ASM" )
write.table(BisSNP_chr20_asm, file = "BisSNP_chr20asm.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#Extract SNPs linked
BisSNP_chr20_snp <- assay(BisSNP_chr20_derASM, "snp.table")
write.table(BisSNP_chr20_snp, file = "BisSNP_chr20_snp.txt" , sep = "\t",
            row.names = TRUE, col.names = NA)
#Extract zASM score
BisSNP_chr20_zASM <- assay(BisSNP_chr20_derASM, "zASM")
write.table(BisSNP_chr20_zASM, file = "BisSNP_chr20_zASM.txt" , sep = "\t",
            row.names = TRUE, col.names = NA)