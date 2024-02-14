#This script runs BWA-meth pipeline on all *.fastq.gz files present within the current working directory.
#The reads are assumed to be named the way they are produced from the sequencing platform.
#Outputs are and can be found in:
#	FastQC (quality control) > ./BWA/FastQC
#	Trimmed reads (using trimgalore) > ./BWA/Trimmed
#	Alignments in SAM format (Default BWA output) > ./BWA/Alignments/Sam
#	Statistics for each alignment > ./BWA/Alignments/Statistics
#	Alignments in BAM format > ./BWA/Alignments/Bam
#	Alignments sorted with samtools > ./BWA/Alignments/Sorted
#	Duplicate marked reads (picard markdup) and flagging metrics > ./BWA/Alignments/MarkDup
#	Methylation levels per CpG > ./BWA/Alignments/MethylDackel/PerCpG
#	Mbias plots > ./BWA/Alignments/MethylDackel/Mbias





echo "This script runs BWA-meth pipeline on all *.fastq.gz files present within the current working directory."

#Make a directory where all of the output will be.
mkdir ./BWA
mkdir ./BWA/Alignments
#FastQC. Quality control step
mkdir ./BWA/FastQC
echo "######################"
echo "                      "
echo "                      "
echo "Starting FastQC!      "
echo "                      "
echo "                      "
echo "######################"

for I in `ls *.fastq.gz` ;
do RawReads=$(basename $I);
FastQCstart=`date +%s`
fastqc -o ./BWA/FastQC ./${RawReads}
FastQCend=`date +%s`
FastQC_RunTime=`expr $FastQCend - $FastQCstart`
echo "FastQC for ${RawReads} took ${FastQC_RunTime} seconds." >> ./BWA/Runtime.txt
done ;

echo "######################"
echo "                      "
echo "                      "
echo "FastQC is done!       "
echo "                      "
echo "                      "
echo "######################"

#Trimming with TrimGalore
mkdir ./BWA/Trimmed
echo "######################"
echo "                      "
echo "                      "
echo "Trimming is starting! "
echo "                      "
echo "                      "
echo "######################"
for I in `ls *_R1_001.fastq.gz` ;
do To_trim=$(basename $I _R1_001.fastq.gz);
Trimstart=`date +%s`
trim_galore --paired -o ./BWA/Trimmed ./${To_trim}_R1_001.fastq.gz ./${To_trim}_R2_001.fastq.gz 2>&1 | tee ./${To_trim}.terminal.txt
Trimend=`date +%s`
Trim_RunTime=`expr $Trimend - $Trimstart`
echo "${To_trim} trimming took ${Trim_RunTime} seconds." >> ./BWA/Partial_Runtime.txt
done ;
cat *terminal.txt > Trimming.terminal_report.txt
rm *terminal.txt
echo "######################"
echo "                      "
echo "                      "
echo "Trimming is done!     "
echo "                      "
echo "                      "
echo "######################"

#Alignment with BWA

mkdir ./BWA/Alignments/Sam
echo "############################"
echo "                            "
echo "                            "
echo "Starting the alignments!    "
echo "                            "
echo "                            "
echo "############################"
for I in `ls ./BWA/Trimmed/*1.fq.gz` ;
do TrimmedReads=$(basename $I R1_001_val_1.fq.gz)
Alignmentstart=`date +%s`
bwameth.py --reference ~/HG38_Combined_with_Ctrl.fa -t 16 ./BWA/Trimmed/${TrimmedReads}R1_001_val_1.fq.gz ./BWA/Trimmed/${TrimmedReads}R2_001_val_2.fq.gz > ./BWA/Alignments/Sam/${TrimmedReads}.sam 2>&1 | tee ./${TrimmedReads}.terminal.txt
Alignmentend=`date +%s`
Alignment_RunTime=`expr $Alignmentend - $Alignmentstart`
echo "${TrimmedReads} alignment took ${Alignment_RunTime} seconds." >> ./BWA/Partial_Runtime.txt
done ;
cat *terminal.txt > Alignment.terminal_report.txt
rm *terminal.txt

echo "######################"
echo "                      "
echo "                      "
echo "Alignments are done!  "
echo "                      "
echo "                      "
echo "######################"

#Get stats for each alignment
cd ./BWA/Alignments
echo "######################"
echo "                      "
echo "                      "
echo "Getting the stats!    "  
echo "                      "
echo "                      "
echo "######################"
for I in `ls ./Sam/*.sam` ;
do AlignmentsSAM=$(basename $I)
Statsstart=`date +%s`
samtools flagstat ./Sam/${AlignmentsSAM} > ./${AlignmentsSAM}.txt 
Statsend=`date +%s`
Stats_RunTime=`expr $Statsend - $Statsstart`
echo "${AlignmentsSAM} statistics took ${Stats_RunTime} seconds." >> ../Partial_Runtime.txt
done ;

mkdir Statistics
mv *.txt Statistics
echo "######################"
echo "                      "
echo "                      "
echo "Stats are here!       "  
echo "                      "
echo "                      "
echo "######################"

#Convert Sam to Bam
mkdir Bam
echo "######################"
echo "                      "
echo "                      "
echo "Converting into BAM!  "
echo "                      "
echo "                      "
echo "######################"
Conversionstart=`date +%s`
samtools view -b -@ 8 ./Sam/${AlignmentsSAM} > ./Bam/${AlignmentsSAM}.bam 2>&1 | tee ./${AlignmentsSAM}.terminal.txt ;
Conversionend=`date +%s`
Conversion_RunTime=`expr $Conversionend - $Conversionstart`
echo "${AlignmentsSAM} conversion to BAM took ${Conversion_RunTime} seconds." >> ../Partial_Runtime.txt
done ;
cat *terminal.txt > Conversion.terminal_report.txt
rm *terminal.txt
echo "######################"
echo "                      "
echo "                      "
echo "Sams are converted!   "
echo "                      "
echo "                      "
echo "######################"

#Sorting
mkdir Sorted

echo "#######################"
echo "                       "
echo "                       "
echo "About to begin sorting!"
echo "                       "
echo "                       "
echo "#######################"

for I in `ls ./Bam/*.bam` ;
do AlignmentsBAM=$(basename $I);
Sortingstart=`date +%s`
samtools sort -@ 8 ./Bam/${AlignmentsBAM} > ./Sorted/sorted.${AlignmentsBAM} 2>&1 | tee ./${AlignmentsBAM}.terminal.txt;
Sortingend=`date +%s`
Sorting_RunTime=`expr $Sortingend - $Sortingstart`
echo "${AlignmentsBAM} sorting took ${Sorting_RunTime} seconds." >> ../Partial_Runtime.txt
done ;
cat *terminal.txt > Sorting.terminal_report.txt
rm *terminal.txt

echo "######################"
echo "                      "
echo "                      "
echo "Sorting is done!"
echo "                      "
echo "                      "
echo "######################"

#Flag Duplicates
mkdir MarkDup
echo "#######################"
echo "                       "
echo "                       "
echo "Marking the duplicates!"
echo "                       "
echo "                       "
echo "#######################"
for I in `ls ./Sorted/*.bam` ;
do SortedReads=$(basename $I);
MarkDupstart=`date +%s`
java -Dpicard.useLegacyParser=false -jar picardcloud.jar MarkDuplicates -I ./Sorted/${SortedReads} -TAG_DUPLICATE_SET_MEMBERS -O ./MarkDup/markdup.${SortedReads} -M ./MarkDup/dupmetrics.${SortedReads} 2>&1 | tee ./${SortedReads}.terminal.txt;
MarkDupend=`date +%s`
MarkDup_RunTime=`expr $MarkDupend - $MarkDupstart`
echo "${SortedReads} marking the duplicates took ${MarkDup_RunTime} seconds." >> ../Partial_Runtime.txt
done ;
cat *terminal.txt > MarkDup.terminal_report.txt
rm *terminal.txt
mkdir ./MarkDup/Metrics
mv ./MarkDup/dupmetrics.* ./MarkDup/Metrics
echo "######################"
echo "                      "
echo "                      "
echo "Duplicates are marked!"
echo "                      "
echo "                      "
echo "######################"

#For MethylDackel, this must be run first in each terminal.
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"
echo "######################"
echo "                      "
echo "                      "
echo "Library is exported"
echo "                      "
echo "                      "
echo "######################"

#Extract Methylation using methyldackel. Some of these can't run simultaneously -mergeContext and -cytosine_report, as well as -counts, -cytosine_report and -methylKit.
mkdir ./MarkDup/MethylDackel
cd ./MarkDup/MethylDackel
echo "################################"
echo "                                "
echo "                                "
echo "The MethylDackel steps begin!   "
echo "                                "
echo "                                "
echo "################################"

echo "The output will be Mbias reports and methylation data per CpG.
	If a different sort of output is desired, either remove hashtags
	from the script next to flags of interest or just include them.
	Some of the possible outputs are cytosine reports and counts.
	For full options, type MethylDackel extract into the terminal 
	and it will provide the manual"


for I in `ls ../*.bam` ;
do MarkedDup=$(basename $I);
MethylDackel extract -cytosine_report -@ 4  ~/HG38_Combined_with_Ctrl.fa ../${MarkedDup};
done ;
#echo "#########################"
#echo "                         "
#echo "                         "
#echo "Cytosine report is done! "
#echo "                         "
#echo "                         "
#echo "#########################"


for I in `ls ../*.bam` ;
do MarkedDup=$(basename $I);
MethylDackel extract -counts -@ 4  ~/HG38_Combined_with_Ctrl.fa ../${MarkedDup};
done ;
#echo "######################"
#echo "                      "
#echo "                      "
#echo "Counts are done!      "
#echo "                      "
#echo "                      "
#echo "######################"



for I in `ls ../*.bam` ;
do MarkedDup=$(basename $I);
Constart=`date +%s`
MethylDackel extract -@ 4 --mergeContext  ~/HG38_Combined_with_Ctrl.fa ../${MarkedDup} 2>&1 | tee ./${MarkedDup}.terminal.txt;
Conend=`date +%s`
Con_RunTime=`expr $Conend - $Constart`
echo "${MarkedDup} took ${Con_RunTime} seconds to merge contexts." >> ./Meth_Runtime.txt
done ;
cat *terminal.txt > mergeContext.terminal_report.txt
rm *terminal.txt
echo "########################"
echo "                        "
echo "                        "
echo "The contexts are merged!"
echo "                        "
echo "                        "
echo "########################"


#Mbias plots
for I in `ls ../*.bam` ;
do MarkedDup=$(basename $I);
Mbiasstart=`date +%s`
MethylDackel mbias  --nOT 6,5,6,5 -@ 4  ~/HG38_Combined_with_Ctrl.fa ../${MarkedDup} ./${MarkedDup}_CpG 2>&1 | tee ./${MarkedDup}.terminal.txt;
Mbiasend=`date +%s`
Mbias_RunTime=`expr $Mbiasend - $Mbiasstart`
echo "${MarkedDup} took ${Mbias_RunTime} to produce Mbias graph." >> ./Meth_Runtime.txt
done ;
cat *terminal.txt > Mbias.terminal_report.txt
rm *terminal.txt
mkdir PerCpG
mkdir Mbias
mv *.bedGraph PerCpG
mv *.svg Mbias
cd ../
cd ../
mv ./MarkDup/MethylDackel ./
cat ../Partial_Runtime.txt ./MethylDackel/Meth_Runtime.txt > BWA_Runtime.txt

echo "#########################"
echo "                         "
echo "                         "
echo "Mbias plots are done!    "
echo "                         "
echo "                         "
echo "#########################"

echo "#########################"
echo "                         "
echo "                         "
echo "The pipeline is finished!"
echo "                         "
echo "                         "
echo "#########################"

echo "Outputs are and can be found in:
	FastQC (quality control) > ./BWA/FastQC
	Trimmed reads (using trimgalore) > ./BWA/Trimmed
	Alignments in SAM format (Default BWA output) > ./BWA/Alignments/Sam
	Statistics for each alignment > ./BWA/Alignments/Statistics
	Alignments in BAM format > ./BWA/Alignments/Bam
	Alignments sorted with samtools > ./BWA/Alignments/Sorted
	Duplicate marked reads (picard markdup) and flagging metrics > ./BWA/Alignments/MarkDup
	Methylation levels per CpG > ./BWA/Alignments/MethylDackel/PerCpG
	Mbias plots > ./BWA/Alignments/MethylDackel/Mbias"







