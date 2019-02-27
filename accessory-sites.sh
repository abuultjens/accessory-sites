
#######################################################################
# accessory-sites.sh
#######################################################################

# sh accessory-sites.sh [INFILE.fa] [PREFIX] [OUTFILE_DATA]

#------------------------------------------------

ALIGNMENT=${1}
#ALIGNMENT=test.fa
#ALIGNMENT=6_sc4.3.6.full.clean.WO-_noref.aln

# display help
if [ "${ALIGNMENT}" == "help" ]; then
    echo ""
    echo "--------------------------------------------------------------------"
    echo "accessory-sites   Andrew Buultjens 2019: buultjensa@gmail.com"
    echo "--------------------------------------------------------------------"
    echo ""
    echo "VERSION:"
    echo "v1.02"
    echo ""
    echo "ABOUT:"
    echo "This program finds invariant accessory sites from a multi FASTA alignment file."
    echo ""
    echo "USAGE:"
    echo "sh accessory-sites.sh [INFILE.fa] [PREFIX] [OUTFILE_DATA]"
    echo ""
    echo "OPTIONS:"
    echo "[OUTFILE_DATA] can be either 'all_sites' or 'only_invariant_accessory_sites'"
    echo "*note that the command options must be in the exact order as specified above as they are treated as positional arguments"    
    echo ""
    echo "EXAMPLE:"
    echo "extracting only invariant accessory sites"
    echo "sh accessory-sites.sh inputfile.aln PREFIX only_invariant_accessory_sites"
    echo ""
    echo "extracting invariant accessory sites and all snp sites"
    echo "sh accessory-sites.sh inputfile.aln PREFIX all_sites"
    echo ""
    
    # crash script and exit
    exit 1
fi

PREFIX=${2}
#PREFIX=OUT

OUTFILE_DATA=${3}
#OUTFILE_DATA=only_invariant_accessory_sites
#OUTFILE_DATA=all_sites

#------------------------------------------------

# generate random prefix for all tmp files
RAND_1=`echo $((1 + RANDOM % 100))`
RAND_2=`echo $((100 + RANDOM % 200))`
RAND_3=`echo $((200 + RANDOM % 300))`
RAND=`echo "${RAND_1}${RAND_2}${RAND_3}"`

#------------------------------------------------
# make files for original alignment

# make files for the original alignment
snp-sites -v ${ALIGNMENT} > ${RAND}_tmp_original.vcf
vcf-to-tab < ${RAND}_tmp_original.vcf > ${RAND}_tmp_original.tab

# get chr col
cat ${RAND}_tmp_original.tab | grep -v "#" | cut -f 1 > ${RAND}_tmp_original_CHR.txt

# get the number of positions in table
original_LC=`cat ${RAND}_tmp_original.tab | grep -v "#" | wc -l | awk '{print $1}'`

# get pos col
cat ${RAND}_tmp_original.tab | grep -v "#" | cut -f 2 > ${RAND}_tmp_original_POS.txt

# make a gff file
seq 1 ${original_LC} | awk '{print "s"}'  > ${RAND}_tmp_original_s.txt
seq 1 ${original_LC} | awk '{print "S"}'  > ${RAND}_tmp_original_S.txt
seq 1 ${original_LC} | awk '{print "."}'  > ${RAND}_tmp_original_dot.txt
seq 1 ${original_LC} | awk '{print "0"}'  > ${RAND}_tmp_original_0.txt
cat ${RAND}_tmp_original.tab | grep -v "#" | cut -f 3- | tr '*' 'N' | tr -d '/' | tr '\t' ',' > ${RAND}_tmp_original_attribute.txt
echo "###gff-version 3" > ${RAND}_tmp_original.gff
paste ${RAND}_tmp_original_CHR.txt ${RAND}_tmp_original_s.txt ${RAND}_tmp_original_S.txt ${RAND}_tmp_original_POS.txt ${RAND}_tmp_original_POS.txt ${RAND}_tmp_original_dot.txt ${RAND}_tmp_original_dot.txt ${RAND}_tmp_original_0.txt ${RAND}_tmp_original_attribute.txt >> ${RAND}_tmp_original.gff

#------------------------------------------------
# make files for fake alignments and compare to original table

# replace N with fake snp
#echo "replacing Ns"
tr 'N' 'A' < ${ALIGNMENT} > ${RAND}_tmp_fake_A.fa &
tr 'N' 'G' < ${ALIGNMENT} > ${RAND}_tmp_fake_G.fa &
tr 'N' 'C' < ${ALIGNMENT} > ${RAND}_tmp_fake_C.fa &
tr 'N' 'T' < ${ALIGNMENT} > ${RAND}_tmp_fake_T.fa &

wait

# run snp-sites
#echo "snp-sites"
snp-sites -v ${RAND}_tmp_fake_A.fa > ${RAND}_tmp_fake_A.vcf &
snp-sites -v ${RAND}_tmp_fake_G.fa > ${RAND}_tmp_fake_G.vcf &
snp-sites -v ${RAND}_tmp_fake_C.fa > ${RAND}_tmp_fake_C.vcf &
snp-sites -v ${RAND}_tmp_fake_T.fa > ${RAND}_tmp_fake_T.vcf &

wait

# run vcf-tools
#echo "vcf tools"
vcf-to-tab < ${RAND}_tmp_fake_A.vcf > ${RAND}_tmp_fake_A.tab &
vcf-to-tab < ${RAND}_tmp_fake_G.vcf > ${RAND}_tmp_fake_G.tab &
vcf-to-tab < ${RAND}_tmp_fake_C.vcf > ${RAND}_tmp_fake_C.tab &
vcf-to-tab < ${RAND}_tmp_fake_T.vcf > ${RAND}_tmp_fake_T.tab &

wait

# if file already exists then delete it to prevent carry over from a previous run
if ls ${RAND}_tmp_found_with_fake.txt 1> /dev/null 2>&1; then
    rm ${RAND}_tmp_found_with_fake.txt
fi

# loop through DNA alphabet
for FAKE in A G C T; do

    # get fake chr col
    cat ${RAND}_tmp_fake_${FAKE}.tab | grep -v "#" | cut -f 1 > ${RAND}_tmp_fake_${FAKE}_CHR.txt
    
    # get the number of positions in table
    fake_LC=`cat ${RAND}_tmp_fake_${FAKE}.tab | grep -v "#" | wc -l | awk '{print $1}'`
    
    # get fake pos col
    cat ${RAND}_tmp_fake_${FAKE}.tab | grep -v "#" | cut -f 2 > ${RAND}_tmp_fake_${FAKE}_POS.txt

    # make a gff file
    seq 1 ${fake_LC} | awk '{print "s"}' > ${RAND}_tmp_fake_${FAKE}_s.txt
    seq 1 ${fake_LC} | awk '{print "S"}' > ${RAND}_tmp_fake_${FAKE}_S.txt
    seq 1 ${fake_LC} | awk '{print "."}' > ${RAND}_tmp_fake_${FAKE}_dot.txt
    seq 1 ${fake_LC} | awk '{print "0"}' > ${RAND}_tmp_fake_${FAKE}_0.txt
    cat ${RAND}_tmp_fake_${FAKE}.tab | grep -v "#" | cut -f 3- | tr '*' 'N' | tr -d '/' | tr '\t' ',' > ${RAND}_tmp_fake_${FAKE}_attribute.txt
    echo "###gff-version 3" > ${RAND}_tmp_fake_${FAKE}.gff
    paste ${RAND}_tmp_fake_${FAKE}_CHR.txt ${RAND}_tmp_fake_${FAKE}_s.txt ${RAND}_tmp_fake_${FAKE}_S.txt ${RAND}_tmp_fake_${FAKE}_POS.txt ${RAND}_tmp_fake_${FAKE}_POS.txt ${RAND}_tmp_fake_${FAKE}_dot.txt ${RAND}_tmp_fake_${FAKE}_dot.txt ${RAND}_tmp_fake_${FAKE}_0.txt ${RAND}_tmp_fake_${FAKE}_attribute.txt >> ${RAND}_tmp_fake_${FAKE}.gff

    # find del variants that are not in the original snp-sites output
    bedtools intersect -v -a ${RAND}_tmp_fake_${FAKE}.gff -b ${RAND}_tmp_original.gff > ${RAND}_tmp_out.gff
    
    # cut off columns of interest and convert fake snps back to Ns
    cat ${RAND}_tmp_out.gff | cut -f 1,4,9 | tr ',' '\t' | tr "${FAKE}" 'N' >> ${RAND}_tmp_found_with_fake.txt

done

# make a non-redundant list of gff hits
sort ${RAND}_tmp_found_with_fake.txt | uniq > ${RAND}_tmp_nr_found_with_fake.txt

# reformat the original data
cat ${RAND}_tmp_original.tab | tr '*' 'N' | tr -d '/' > ${RAND}_tmp_original_tmp.tab

#-----------------------------------------------------
# select which data to include in outfile

# output only invariant accessory sites
if [ "${OUTFILE_DATA}" == "only_invariant_accessory_sites" ]; then
    cat ${RAND}_tmp_original_tmp.tab | grep "#" > ${RAND}_tmp_HEAD.txt
    cat ${RAND}_tmp_HEAD.txt ${RAND}_tmp_nr_found_with_fake.txt | cut -f 1,2,4- > ${RAND}_tmp_${PREFIX}.tab
fi

# output snps and invariant accessory sites
if [ "${OUTFILE_DATA}" == "all_sites" ]; then
    cat ${RAND}_tmp_original_tmp.tab ${RAND}_tmp_nr_found_with_fake.txt | cut -f 1,2,4- > ${RAND}_tmp_${PREFIX}.tab
fi

#-----------------------------------------------------
# make outfiles

# output table
cp ${RAND}_tmp_${PREFIX}.tab ${PREFIX}.tab

# output fasta
NAMES_COL=`cat ${RAND}_tmp_${PREFIX}.tab | datamash transpose -H | grep -v "CHROM" | grep -v "POS" | cut -f 1 > ${RAND}_tmp_NAMES_COL.txt`
SEQ_COL=`cat ${RAND}_tmp_${PREFIX}.tab | datamash transpose -H | grep -v "CHROM" | grep -v "POS" | cut -f 2- | tr -d '\t' > ${RAND}_tmp_SEQ_COL.txt`
paste ${RAND}_tmp_NAMES_COL.txt ${RAND}_tmp_SEQ_COL.txt | awk '$0=">"$0' | tr '\t' '\n' > ${PREFIX}.fa


#-----------------------------------------------------

# remove all tmp files
rm ${RAND}_tmp_*

#-----------------------------------------------------

