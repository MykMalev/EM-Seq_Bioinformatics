#Reformat DAMEfinder calc_derivedasm output's snp.table into a bed format.
#Two outputs are produced: "SNPs.LinkedToCpG.bedgraph" is a list of SNPs (chr and position) that were considered in the ASM analysis. "snps.and.cpgs.bedgraph" is a list of SNPs (chr and position) together with CpGs linked to them.
makesnp()
{
sed -i '1d' $1
awk -v OFS='\t' '{ print $1 }' $1 > $1.tmp
awk -v OFS='\t' '{ print $2 }' $1 > $1.tmp2
cut -d '.' -f 2 <$1.tmp2 >$1.tmp3
cut -d '.' -f 2 <$1.tmp3 >$1.tmp4
awk -F'\t' -vOFS='\t' '{ $1 = "chr" $1 }1' $1.tmp > $1.tmp5
cut -d '.' -f 1 <$1.tmp5 > $1.tmp12
awk -F'\t' '{gsub(/"/, "", $1); print $1}' $1.tmp12 > $1.tmp8
awk -F'\t' '{gsub(/"/, "", $1); print $1}' $1.tmp4 > $1.tmp6
awk -F'\t' '{gsub(/"/, "", $1); print $1}' $1.tmp3 > $1.tmp7
awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],$0}' $1.tmp8 $1.tmp8 > $1.tmp9
awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],$0}' $1.tmp9 $1.tmp7 > $1.tmp10
awk -v OFS='\t' ' {$4=$2+1}{ print $1, $2, $3, $4 }' $1.tmp10 > $1.tmp11
awk -v OFS='\t' '{ print $1, $2 , $4, $3 }' $1.tmp11 > $1.snps.and.cpgs.bedgraph
awk -v OFS='\t' '{ print $1, $4 }' $1.snps.linked.to.cpgs > $1.tmp14
awk -v OFS='\t' '{$3=$2+1} {$4=1} { print $1, $2, $3, $4 }' $1.tmp14 > $1.tmp13
awk '!a[$0]++' $1.tmp13 > $1.SNPs.LinkedToCpG.bedgraph
rm $1.tmp
rm $1.tmp2
rm $1.tmp3
rm $1.tmp4
rm $1.tmp5
rm $1.tmp6
rm $1.tmp7
rm $1.tmp8
rm $1.tmp9
rm $1.tmp10
rm $1.tmp11
rm $1.tmp12
rm $1.tmp13
rm $1.tmp14 
}

### For a specific input
#echo "makesnp"
#makesnp $1

### For all text files in the folder
echo "makesnp"
for I in `ls *snp*` ;
do 
makesnp $I
done ;

