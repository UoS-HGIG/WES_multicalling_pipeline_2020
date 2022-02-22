vcftools --vcf $1 --extract-FORMAT-info GQ --out temp0

awk '{a=0; for(i = 3; i<= NF; i++) {a=a+$i}; {print a/NF}}' temp0.GQ.FORMAT > temp_meanGQ.out

grep "^##" $1 >temp_header
grep -v "^##" $1 >temp_body

sed -i "1s|0|99.0|" temp_meanGQ.out

paste temp_meanGQ.out temp_body > temp1

awk '{ if($1>35) {print $0}}' temp1 >temp2

cut -f 2- temp2 >temp3

cat temp_header temp3 > $2
