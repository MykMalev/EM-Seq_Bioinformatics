#####TUPLE based

#DAMEFinder
#First, a methtuple file must be created by using MethTuple. --sc refers to "strand collapse", -m 2 is tuple size. These settings were selected according to DAMEFinder requirements.
methtuple --sc -m 2 -o ${BAM}.tuple --gzip $BAM

#In R
#Load libraries required
library(DAMEfinder)
library(GenomicRanges)
library(SummarizedExperiment)

#Set functions for importing files
DATA_PATH_DIR <- system.file('extdata', package = 'DAMEfinder')
get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)
ref <- get_data_path('Combined_Genome.fa') #Reference genome

#Import the methtuple output.
TupleFile <- get_data_path("TupleFile.tsv.gz")

#read_tuples is a DAMEfinder command that imports the tuple file and produces a list of tibbles.
ReadTuple <- read_tuples(TupleFile, c("SampleName"), minCoverage = 6, maxGap = 20, verbose = TRUE)

#calc_asm calculates ASM values from the read_tuples generated object
TupleCalcASM <- calc_asm(
  ReadTuple,
  beta = 0.5,
  a = 0.2,
  coverage = 5,
  verbose = TRUE
)

#In order to retrieve asm values from the object, "assay" command is used. In order to view all of the information available, "assays" can be used.
TupleASM <- assay(TupleCalcASM, "asm")

#Finally, the tuple based asm scores can be exported as a tab separated table.
write.table(TupleASM, file = "TupleASM", sep = "\t",
            row.names = TRUE, col.names = NA)
            
#Outside R
#The generated table is poorly formated, in order to convert it into a bedgraph that can be viewed in IGV or formated, I wrote this. I usually use it in a while loop.
ls *ASM > files

tupletobed()
{
sed -i '1d' $1 ;	# remove header. sed removes the first line of the original file, so make sure to have a backup or remove this line if the code were to be used on the same file again
awk -F'\t' '{gsub(/"/, "", $1); print $1, $2}' $1 > $1.tmp ; # remove quotation marks
awk -v OFS='\t' '{ print $1 }' $1.tmp > $1.tmp2 #extract first col
awk -F'.' '{ print $1"\t"$2"\t"$3 }' $1.tmp2 > $1.tmp3 # tab  separate it
awk -v OFS='\t' '{print $2}' $1 > $1.tmp4 #take out col2
tail -n +1 $1.tmp4 > $1.tmp5 #remove the first line which is just empty
awk -F'\t' 'NR==FNR{a[NR]=$0;next}{print a[FNR],$0}' $1.tmp3 $1.tmp5 > $1.tuple.bedgraph #combine both files
rm $1.tmp
rm $1.tmp2
rm $1.tmp3
rm $1.tmp4
rm $1.tmp5
}

while read line
do
tupletobed $line
done < files

#DNMTools
#First step is running "states". The output is in epiread format which will be used as an input for later steps.
/home/student/Documents/Software/dnmtools-1.4.2/dnmtools states -z  -o ~/States -t 4 -c ~/reference_genome  ~/mapped_reads

#Allelic outputs allele specifically methylated sites
/home/student/Documents/Software/dnmtools-1.4.2/dnmtools allelic -v -o ./ASM.Sites -c ~/reference_genome ~/States

#AMR outputs allele specifically methylated regions
/home/student/Documents/Software/dnmtools-1.4.2/dnmtools amrfinder -o ./ASM.Regions -c ~/reference_genome ~/States

#AMR may fail and outpute the error "Cannot allocate memory". In that case, the States file can be separated in batches of chromosomes using grep. Once the AMR analysis is finished, the output can be merged using cat.

head chr
chr1
chr2
chr3
chr4

grep -w --file=chr States > chr1to4.States


#####SNP-based

#CGmapTools
#Allele specifically methylated sites
cgmaptools asm -o ~/CGmapTools.ASS -m ass  -r ~/reference_genome -b ~/mapped_reads -l ~/SNP_File

#Allele specifically methylated regions
cgmaptools asm -o ~/CGmapTools.ASM  -r ~/reference_genome -b ~/mapped_reads -l ~/SNP_File

#DAMEfinder
#In R
#Load libraries
library(DAMEfinder)
library(GenomicRanges)
library(SummarizedExperiment)
#File retrieval functions
DATA_PATH_DIR <- system.file('extdata', package = 'DAMEfinder')
get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)
ref <- get_data_path('Combined_Genome.fa')

#Import mapped reads and SNPs
samplebam <- get_data_path("sample.bam")
samplevcf <- get_data_path("sample.vcf")
sampleNames <- "SampleName"
#extract_bams extracts reads from the bam file according to alleles and generates a list in GRangesList format
extractedsample <- extract_bams(samplebam, samplevcf, sampleNames, coverage = 6, ref, cores = 32)

#calc_derivedasm calculates ASM values and generates a RangedSummarizedExperiment object.
SampleDerivedASM <- calc_derivedasm(extractedsample, cores = 32, verbose = TRUE)

#Extract ASM scores and write the table.
SampleASM <- assay(SampleDerivedASM, "der.ASM" )
write.table(SampleASM, file = "SampleASM", sep = "\t",
            row.names = TRUE, col.names = NA)
            
#Outside R
#The ASM table generated is formated poorly, thus it is converted to a bedgraph with the following code. 
asmtobed()
{
sed -i '1d' $1 ; #removes header
awk -F'\t' '{gsub(/"/, "", $1); print $1, $2}' $1 > $1.tmp ; #Removes quotation marks from the first column
cut -d '.' -f 2 <$1.tmp > $1.postest ; #Removes the dots
awk -v OFS='\t' '{ print $1 }' $1.postest > $1.positions ; #Extracts only the first column
awk -v OFS='\t' ' {$2=$1+1}{ print $1, $2 }' $1.positions > $1.startend ; # Creates start-end coordinates (Adds 1 to a positon)
awk -F'\t' -vOFS='\t' '{ $1 = "chr" $1 }1' $1.tmp > $1.tmp2 ; #Adds "chr" next to position
awk '{sub(/\..*$/, "", $1)} 1' $1.tmp2 > $1.tmp3 ; #Removes the coordinates
awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],$0}' $1.tmp3 $1.startend > $1.full ; # Combines the file with chr and score together with start-end coordinates
awk -v OFS='\t' '{ print $1, $3 , $4, $2 }' $1.full > $1.derived.asm.bedgraph ; #Rearranges the columns like in a typical bedgraph format
rm $1.tmp
rm $1.postest
rm $1.positions
rm $1.startend
rm $1.tmp2
rm $1.tmp3
rm $1.full
}

asmtobed() $1
