#Determine EM-Seq conversion rates by analysing puc19 and lambda methylation status from CpG methylation coverage bedgraphs

grep -w "M77789.2" ~/CpG.bedGraph > puc19


grep -w "J02459.1" ~/CpG.bedGraph > lambda


awk '{ total += $4 } END { print total/NR }' puc19 > puc19.scoreavg 

awk '{ total += $4 } END { print total/NR }' lambda > lambda.scoreavg 
