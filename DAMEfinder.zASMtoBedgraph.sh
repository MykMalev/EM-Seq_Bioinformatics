#Reformat DAMEfinder calc_derivedasm output's zASM scores into a bedgraph format.

zasmtobed ()
{
sed -i '1d' $1 ;
awk -F'\t' '{gsub(/"/, "", $1); print $1, $2}' $1 > $1.tmp ;
cut -d '.' -f 2 <$1.tmp > $1.tmp2 ;
awk -v OFS='\t' '{ print $1 }' $1.tmp2 > $1.tmp3 ;
awk -v OFS='\t' ' {$2=$1+1}{ print $1, $2 }' $1.tmp3 > $1.tmp4 ;
awk -F'\t' -vOFS='\t' '{ $1 = "chr" $1 }1' $1.tmp > $1.tmp5 ;
awk '{sub(/\..*$/, "", $1)} 1' $1.tmp2 > $1.tmp6 ;
awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],$0}' $1.tmp3 $1.tmp4 > $1.tmp7 ;
awk -v OFS='\t' '{ print $1, $3 , $4, $2 }' $1.tmp7 > $1.tmp8
sed '/NA/d' $1.zASMtmp > $1.zASM.bedgraph ;
rm $1.tmp
rm $1.tmp2
rm $1.tmp3 
rm $1.tmp4
rm $1.tmp5
rm $1.tmp6
rm $1.tmp7
rm $1.tmp8
}

### For every text file in folder
echo "zASMtoBed"
for I in `ls *zASM*` ;
do 
zasmtobed $I
done ;

## For a specific input
#zasmtobed $1
