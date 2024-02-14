#This script runs the Bismark alignment pipeline for all *.fastq.gz files within the current working directory.
#The reads are assumed to be named the way they are produced from the sequencing platform.
#Outputs are and can be found in:
#	FastQC (quality control) > ./Bismark/FastQC
#	Trimmed reads (using trimgalore) > ./Bismark/Trimmed
#	Alignments > ./Bismark/Alignments
#	Deduplicated alignments and deduplications reports > ./Bismark/Deduplicated
#	Methylation data (Splitting Reports, Mbias plots, coverage reports, bedrgaphs, cytosine context summaries, CpG contexts) > ./Bismark/MethylationData

echo "This script runs the Bismark alignment pipeline for all *.fastq.gz files within the current working directory"

#Make a directory for the output.

mkdir ./Bismark


#FastQC. Quality control step.

mkdir ./Bismark/FastQC
for I in `ls *.fastq.gz` ;
do RawReads=$(basename $I);
FastQCstart=`date +%s`
fastqc -o ./Bismark/FastQC ./${RawReads}
FastQCend=`date +%s`
FastQC_runtime=`expr $FastQCend - $FastQCstart`
echo "FastQC for ${RawReads} took ${FastQC_runtime} seconds." >> ./Bismark/Partial_Runtime.txt
done ;


echo "###############"
echo "FastQC is done!"
echo "###############"

#Trimming with TrimGalore.

mkdir ./Bismark/Trimmed
for I in `ls *_R1_001.fastq.gz` ;
do To_trim=$(basename $I _R1_001.fastq.gz);
Trimstart=`date +%s`
trim_galore --paired -o ./Bismark/Trimmed ./${To_trim}_R1_001.fastq.gz ./${To_trim}_R2_001.fastq.gz 2>&1 | tee ./${To_trim}.terminal.txt
Trimend=`date +%s`
Trim_runtime=`expr $Trimend - $Trimstart`
echo "Trimming for ${To_trim} took ${Trim_runtime} seconds." >> ./Bismark/Partial_Runtime.txt
done ;
cat *terminal.txt > TrimGalore.terminal_report.txt
rm *terminal.txt

echo "#################"
echo "Trimming is done!"
echo "#################"

#Alignment with Bismark.

mkdir ./Bismark/Alignments
mkdir ./Bismark/Alignments/UnmappedReads
mkdir ./Bismark/Alignments/AmbiguousReads
mkdir ./Bismark/Alignments/AlignmentReport

for I in `ls ./Bismark/Trimmed/*1.fastq.gz` ;
do TrimmedReads=$(basename $I R1_001_val_1.fastq.gz)
Alignmentstart=`date +%s`
bismark --path_to_bowtie ~/bowtie2-2.4.2-linux-x86_64 --samtools_path ~/samtools --genome ~/HG38_with_Ctrl/ --dovetail --parallel 6 --unmapped --ambiguous -o ./Bismark/Alignments --gzip -1 ./Bismark/Trimmed/${TrimmedReads}R1_001_val_1.fastq.gz -2 ./Bismark/Trimmed/${TrimmedReads}R2_001_val_2.fastq.gz 2>&1 | tee ./${TrimmedReads}.terminal.txt
Alignmentend=`date +%s`
Alignment_runtime=`expr $Alignmentend - $Alignmentstart`
echo "Alignment for ${TrimmedReads} took ${Alignment_runtime} seconds." >> ./Bismark/Partial_Runtime.txt
done ;
cat *terminal.txt > Alignment.terminal_report.txt
rm *terminal.txt

mv ./Bismark/Alignments/*unmapped* ./Bismark/Alignments/UnmappedReads
mv ./Bismark/Alignments/*ambiguous* ./Bismark/Alignments/AmbiguousReads
mv ./Bismark/Alignments/*.txt ./Bismark/Alignments/AlignmentReport

echo "####################"
echo "Alignments are done!"
echo "####################"

#Removing the duplicate reads.

cd ./Bismark/Alignments
mkdir Deduplicated
mkdir Deduplicated/DedupReports
for I in `ls *.bam` ;
do Alignment=$(basename $I)
Dedupstart=`date +%s`
~/deduplicate_bismark --output_dir ./Deduplicated ./${Alignment} 2>&1 | tee ./${Alignment}.terminal.txt
Dedupend=`date +%s`
Dedup_runtime=`expr $Dedupend - $Dedupstart`
echo "${Alignment} deduplication took ${Dedup_runtime} seconds." >> ../Partial_Runtime.txt
done ;
cat *terminal.txt > Deduplication.terminal_report.txt
rm *terminal.txt

mv ./Deduplicated/*.txt ./Deduplicated/DedupReports
echo "######################"
echo "Deduplication is done!"
echo "######################"
#Extract methylation data.

mkdir MethylationData

for I in `ls ./Deduplicated/*.bam` ;
do Deduplicates=$(basename $I)
Methstart=`date +%s`
~/bismark_methylation_extractor --comprehensive --merge_non_CpG --bedgraph --cytosine_report --multicore 8 --ignore_r2 5 --ignore 5 --ignore_3prime_r2 5 --ignore_3prime 5 --samtools_path ~/samtools --genome_folder ~/HG38_with_Ctrl/ -o ./MethylationData --gzip ./Deduplicated/${Deduplicates} 2>&1 | tee ./${Deduplicates}.terminal.txt
Methend=`date +%s`
Meth_runtime=`expr $Methend - $Methstart`
echo "${Deduplicates} methylation extraction took ${Meth_runtime} seconds." >> ../Partial_Runtime.txt
done ;
cat *terminal.txt > MethylationData.terminal_report.txt
rm *terminal.txt



echo "###############################"
echo "Methylation extraction is done!"
echo "###############################"
#Merge both strands.
cd ./MethylationData
mkdir CombinedCpG

for I in `ls ./*bismark.cov.gz` ;
do CoverageFile=$(basename $I)
Mergestart=`date +%s`
~/coverage2cytosine --genome_folder ~/HG38_with_Ctrl/ --merge_CpG --dir ./CombinedCpG -o ${CoverageFile} --gzip ./${CoverageFile} 2>&1 | tee ./${CoverageFile}.terminal.txt
Mergeend=`date +%s`
Merge_runtime=`expr $Mergeend - $Mergestart`
echo "${CoverageFile} merging both strands took ${Merge_runtime} seconds." >> ./Merge_Runtime.txt
done ;
cat *terminal.txt > CoverageFile.terminal_report.txt
rm *terminal.txt



echo "#########################"
echo "Both strands are merged!"
echo "#########################"
#Coverage report.

mkdir CoverageReport

for I in `ls ../Deduplicated/*.bam` ;
do Deduplicates=$(basename $I)
Covstart=`date +%s`
~/bam2nuc --genome_folder ~/HG38_with_Ctrl/ --samtools_path ~/samtools --dir ./CoverageReport ../Deduplicated/${Deduplicates} 2>&1 | tee ./${Deduplicates}.terminal.txt
Covend=`date +%s`
Cov_runtime=`expr $Covend - $Covstart`
echo "${Deduplicates} generating coverage file took ${Cov_runtime} seconds." >> ./Merge_Runtime.txt
done ;
cat *terminal.txt > CoverageReport.terminal_report.txt
rm *terminal.txt




mkdir bedGraphs
mkdir CpG_Context
mkdir MbiasTextFiles
mkdir SplittingReports
mkdir CytosineContext
mkdir MbiasPlots

mv ./Merge_Runtime.txt ../
mv *bedGraph* bedGraphs
mv *CpG_context_* CpG_Context
mv *M-bias.txt MbiasTextFiles
mv *splitting_* SplittingReports
mv *cytosine_* CytosineContext
mv *.png MbiasPlots
mv *CpG_report.txt.gz CombinedCpG
cd ../
cat ../Partial_Runtime.txt ./Merge_Runtime.txt > Bismark_Runtime.txt

echo "###########################"
echo "Coverage reports are done!"
echo "###########################"

echo "#################################"
echo "The Bismark pipeline is finished!"
echo "#################################"

echo "Outputs are and can be found in:
	FastQC (quality control) > ./Bismark/FastQC
	Trimmed reads (using trimgalore) > ./Bismark/Trimmed
	Alignments > ./Bismark/Alignments
	Deduplicated alignments and deduplications reports > ./Bismark/Deduplicated
	Methylation data (Splitting Reports, Mbias plots, coverage reports, bedrgaphs, cytosine context summaries, CpG contexts) > ./Bismark/MethylationData"




