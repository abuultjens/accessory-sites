
#######################################################################
# accessory-sites.sh
#######################################################################

# sh accessory-sites.sh [aln.fa] [prefix]

ALIGNMENT=${1}
PREFIX=${2}

# make files for the original alignment
snp-sites -v ${ALIGNMENT} > original.vcf
vcf-to-tab < original.vcf > original.tab

# get chr col
cat original.tab | grep -v "#" | cut -f 1 > original_CHR.txt

# get pos col
cat original.tab | grep -v "#" | cut -f 2 > original_POS.txt

# make a gff file
tr '1' 's' < original_CHR.txt > original_s.txt
tr '1' 'S' < original_CHR.txt > original_S.txt
tr '1' '.' < original_CHR.txt > original_dot.txt
tr '1' '0' < original_CHR.txt > original_0.txt
cat original.tab | grep -v "#" | cut -f 3- | tr '*' 'N' | tr -d '/' | tr '\t' ',' > original_attribute.txt
echo "###gff-version 3" > original.gff
paste original_CHR.txt original_s.txt original_S.txt original_POS.txt original_POS.txt original_dot.txt original_dot.txt original_0.txt original_attribute.txt >> original.gff

# replace N with fake snp
#echo "replacing Ns"
tr 'N' 'A' < ${ALIGNMENT} > fake_A.fa &
tr 'N' 'G' < ${ALIGNMENT} > fake_G.fa &
tr 'N' 'C' < ${ALIGNMENT} > fake_C.fa &
tr 'N' 'T' < ${ALIGNMENT} > fake_T.fa &

wait

# run snp-sites
#echo "snp-sites"
snp-sites -v fake_A.fa > fake_A.vcf &
snp-sites -v fake_G.fa > fake_G.vcf &
snp-sites -v fake_C.fa > fake_C.vcf &
snp-sites -v fake_T.fa > fake_T.vcf &

wait

# run vcf-tools
#echo "vcf tools"
vcf-to-tab < fake_A.vcf > fake_A.tab &
vcf-to-tab < fake_G.vcf > fake_G.tab &
vcf-to-tab < fake_C.vcf > fake_C.tab &
vcf-to-tab < fake_T.vcf > fake_T.tab &

wait



# if file already exists then delete it to prevent carry over from a previous run
if ls found_with_fake.txt 1> /dev/null 2>&1; then
    rm found_with_fake.txt
fi

# loop through DNA alphabet
for FAKE in A G C T; do

    # get fake chr col
    cat fake_${FAKE}.tab | grep -v "#" | cut -f 1 > fake_${FAKE}_CHR.txt
    
    # get fake pos col
    cat fake_${FAKE}.tab | grep -v "#" | cut -f 2 > fake_${FAKE}_POS.txt

    # make a gff file
    tr '1' 's' < fake_${FAKE}_CHR.txt > fake_${FAKE}_s.txt
    tr '1' 'S' < fake_${FAKE}_CHR.txt > fake_${FAKE}_S.txt
    tr '1' '.' < fake_${FAKE}_CHR.txt > fake_${FAKE}_dot.txt
    tr '1' '0' < fake_${FAKE}_CHR.txt > fake_${FAKE}_0.txt
    cat fake_${FAKE}.tab | grep -v "#" | cut -f 3- | tr '*' 'N' | tr -d '/' | tr '\t' ',' > fake_${FAKE}_attribute.txt
    echo "###gff-version 3" > fake_${FAKE}.gff
    paste fake_${FAKE}_CHR.txt fake_${FAKE}_s.txt fake_${FAKE}_S.txt fake_${FAKE}_POS.txt fake_${FAKE}_POS.txt fake_${FAKE}_dot.txt fake_${FAKE}_dot.txt fake_${FAKE}_0.txt fake_${FAKE}_attribute.txt >> fake_${FAKE}.gff

    # find del variants that are not in the original snp-sites output
    bedtools intersect -v -a fake_${FAKE}.gff -b original.gff > out.gff
    
    # cut off columns of interest and convert fake snps back to Ns
    cat out.gff | cut -f 1,4,9 | tr ',' '\t' | tr "${FAKE}" 'N' >> found_with_fake.txt

done

# make a non-redundant list of gff hits
sort found_with_fake.txt | uniq > nr_found_with_fake.txt

# reformat the original fasta data
cat original.tab | tr '*' 'N' | tr -d '/' > original_tmp.tab

# cat the snps and dels to outfile
cat original_tmp.tab nr_found_with_fake.txt | cut -f 1,2,4- > ${PREFIX}.tab

N_SNPS=`cat original_tmp.tab | grep -v "#" | wc -l`
N_DELS=`cat nr_found_with_fake.txt | wc -l`
#echo ""
#echo "--------------------------------"
#echo "summary:"
#echo "found ${N_SNPS} SNPs"
#echo "found ${N_DELS} deletions"
#echo "--------------------------------"

datamash transpose -H < ${PREFIX}.tab | grep -v "#" | grep -v "POS" > ${PREFIX}.tr.tab

#echo ""
#echo "--------------------------------"
#echo "outfile:"
cat ${PREFIX}.tab | datamash transpose -H | grep -v "CHROM"