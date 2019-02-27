#!/bin/bash

#######################################################################
# accessory-sites.sh
#######################################################################

# sh accessory-sites.sh [INFILE.fa] [PREFIX] [OUTFILE_DATA]

#------------------------------------------------

ALIGNMENT=${1}
#ALIGNMENT=test.fa
#ALIGNMENT=test_no-var.fa
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

echo "Options:" > log.txt
echo "ALIGNMENT: ${ALIGNMENT}" >> log.txt
echo "PREFIX: ${PREFIX}" >> log.txt
echo "OUTFILE_DATA: ${OUTFILE_DATA}" >> log.txt

#------------------------------------------------

# check that infile exists
if [ -e "${ALIGNMENT}" ] 
then
    echo "infile exitst" >> log.txt
else
    # print error
    echo "ERROR: cannot find ${ALIGNMENT}"  
    echo "infile does not exitst" >> log.txt
    # crash script and exit
    exit 1
fi

#------------------------------------------------

# generate random prefix for all tmp files
RAND_1=`echo $((1 + RANDOM % 100))`
RAND_2=`echo $((100 + RANDOM % 200))`
RAND_3=`echo $((200 + RANDOM % 300))`
RAND=`echo "${RAND_1}${RAND_2}${RAND_3}"`
echo "random prefix: ${RAND}" >> log.txt

#------------------------------------------------

# make files for the original alignment
snp-sites -v ${ALIGNMENT} > ${RAND}_tmp_original.vcf

# check if snps were found in original alignment
N_SNPS=`cat ${RAND}_tmp_original.vcf | grep -v "#"| wc -l | awk '{print $1}'`
    
####### IF NO CORE SNPS
if [ "$N_SNPS" == "0" ]; then
    # print error
    echo "No core snps were found"  

####### IF YES CORE SNPS
else
    vcf-to-tab < ${RAND}_tmp_original.vcf > ${RAND}_tmp_original.tab
    TMP=`wc -l ${RAND}_tmp_original.tab | awk '{print $1}'`
    echo "made ${RAND}_tmp_original.tab, ${TMP} lines" >> log.txt
    
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

fi

#------------------------------------------------

# make files for fake alignments and compare to original table

# replace N with fake snp
tr 'N' 'A' < ${ALIGNMENT} > ${RAND}_tmp_fake_A.fa &
tr 'N' 'G' < ${ALIGNMENT} > ${RAND}_tmp_fake_G.fa &
tr 'N' 'C' < ${ALIGNMENT} > ${RAND}_tmp_fake_C.fa &
tr 'N' 'T' < ${ALIGNMENT} > ${RAND}_tmp_fake_T.fa &
echo "made fake alignments" >> log.txt


wait

# run snp-sites
snp-sites -v ${RAND}_tmp_fake_A.fa > ${RAND}_tmp_fake_A.vcf &
snp-sites -v ${RAND}_tmp_fake_G.fa > ${RAND}_tmp_fake_G.vcf &
snp-sites -v ${RAND}_tmp_fake_C.fa > ${RAND}_tmp_fake_C.vcf &
snp-sites -v ${RAND}_tmp_fake_T.fa > ${RAND}_tmp_fake_T.vcf &

wait

# check if snps were found
WC_A=`wc -l ${RAND}_tmp_fake_A.vcf | awk '{print $1}'`
echo "made ${RAND}_tmp_fake_A.vcf, ${WC_A} lines" >> log.txt
WC_G=`wc -l ${RAND}_tmp_fake_G.vcf | awk '{print $1}'`
echo "made ${RAND}_tmp_fake_G.vcf, ${WC_G} lines" >> log.txt
WC_C=`wc -l ${RAND}_tmp_fake_C.vcf | awk '{print $1}'`
echo "made ${RAND}_tmp_fake_C.vcf, ${WC_C} lines" >> log.txt
WC_T=`wc -l ${RAND}_tmp_fake_T.vcf | awk '{print $1}'`
echo "made ${RAND}_tmp_fake_T.vcf, ${WC_T} lines" >> log.txt

wait

# run vcf-tools if snps were found
if [ "$WC_A" -gt "0" ]; then
    vcf-to-tab < ${RAND}_tmp_fake_A.vcf > ${RAND}_tmp_fake_A.tab &
    cat ${RAND}_tmp_fake_A.tab | grep "#" > ${RAND}_tmp_HEAD.txt
fi
if [ "$WC_G" -gt "0" ]; then
    vcf-to-tab < ${RAND}_tmp_fake_G.vcf > ${RAND}_tmp_fake_G.tab &
    cat ${RAND}_tmp_fake_G.tab | grep "#" > ${RAND}_tmp_HEAD.txt
fi
if [ "$WC_C" -gt "0" ]; then
    vcf-to-tab < ${RAND}_tmp_fake_C.vcf > ${RAND}_tmp_fake_C.tab &
    cat ${RAND}_tmp_fake_C.tab | grep "#" > ${RAND}_tmp_HEAD.txt
fi
if [ "$WC_T" -gt "0" ]; then
    vcf-to-tab < ${RAND}_tmp_fake_T.vcf > ${RAND}_tmp_fake_T.tab & 
fi

wait

if [ ! -e "${RAND}_tmp_fake_A.tab" ]; then
    if [ ! -e "${RAND}_tmp_fake_G.tab" ]; then
        if [ ! -e "${RAND}_tmp_fake_C.tab" ]; then
            if [ ! -e "${RAND}_tmp_fake_T.tab" ]; then
                echo "*** No invariant acessory sites were found ***"
                exit 1
            fi
        fi
    fi
fi

# if file already exists then delete it to prevent carry over from a previous run
if ls ${RAND}_tmp_found_with_fake.txt 1> /dev/null 2>&1; then
    rm ${RAND}_tmp_found_with_fake.txt
fi

# loop through DNA alphabet
for FAKE in A G C T; do

    ####### IF INVARIANT ACCESSORY SITES FOUND
    if [ -e "${RAND}_tmp_fake_${FAKE}.tab" ]; then

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

        ####### IF YES CORE SNPS
        if [ -e "${RAND}_tmp_original.tab" ]; then
            # find del variants that are not in the original snp-sites output
            bedtools intersect -v -a ${RAND}_tmp_fake_${FAKE}.gff -b ${RAND}_tmp_original.gff > ${RAND}_tmp_out.gff
        
        ####### IF NO CORE SNPS
        else
            # if no snps were found skip bedtools intercept and treat all findings as novel
            mv ${RAND}_tmp_fake_${FAKE}.gff ${RAND}_tmp_out.gff
        fi
    
        # cut off columns of interest and convert fake snps back to Ns
        cat ${RAND}_tmp_out.gff | cut -f 1,4,9 | tr ',' '\t' | tr "${FAKE}" 'N' >> ${RAND}_tmp_found_with_fake.txt
    
    fi

done

# make a non-redundant list of gff hits
sort ${RAND}_tmp_found_with_fake.txt | uniq | grep -v "###gff-version" > ${RAND}_tmp_nr_found_with_fake.txt

# count how many invariant accessory sites were found
N_IAC=`cat ${RAND}_tmp_nr_found_with_fake.txt| wc -l | awk '{print $1}'`

# make header for output
cat ${RAND}_tmp_fake_?.tab | grep "#" | head -1 > ${RAND}_tmp_HEAD.txt

####### IF YES CORE SNPS
if [ -e "${RAND}_tmp_original.tab" ]; then
    # reformat the original data
    cat ${RAND}_tmp_original.tab | tr '*' 'N' | tr -d '/' > ${RAND}_tmp_original_tmp.tab
fi

#-----------------------------------------------------

# select which data to include in outfile

####### IF YES CORE SNPS
if [ -e "${RAND}_tmp_original.tab" ]; then
    # output only invariant accessory sites
    if [ "${OUTFILE_DATA}" == "only_invariant_accessory_sites" ]; then
        cat ${RAND}_tmp_HEAD.txt ${RAND}_tmp_nr_found_with_fake.txt | cut -f 1,2,4- > ${RAND}_tmp_${PREFIX}.tab
    fi

    # output snps and invariant accessory sites
    if [ "${OUTFILE_DATA}" == "all_sites" ]; then
        cat ${RAND}_tmp_original_tmp.tab ${RAND}_tmp_nr_found_with_fake.txt | cut -f 1,2,4- > ${RAND}_tmp_${PREFIX}.tab
    fi
####### IF NO CORE SNPS
else
    cat ${RAND}_tmp_HEAD.txt ${RAND}_tmp_nr_found_with_fake.txt | cut -f 1,2,4- > ${RAND}_tmp_${PREFIX}.tab
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

echo "Found ${N_SNPS} snp sites"
echo "Found ${N_IAC} invariant accessory sites"


