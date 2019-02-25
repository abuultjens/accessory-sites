# accessory-sites
A script to extract accessory sites from whole genome alignments

# Usage
$ sh accessory-sites.sh [INFILE.fa] [PREFIX] [OUTFILE_FORMAT] [OUTFILE_DATA]   

# Options
[OUTFILE_FORMAT] can be either 'table' or 'fasta'   
[OUTFILE_DATA] can be either 'all_sites' or 'only_invariant_accessory_sites'

# How it works
accessory-sites tricks the incredibly fast snp-sites into thinking that invariant accessory positions actually contain snps. This is useful as snp-sites then rapidly extracts the fake snp sites which are then converted back to being invariant accessory sites and concatinated to the original variant core and accessory sites extracted by snp-sites.

# Dependencies
vcf-tools     
bedtools     
snp-sites    

# Example

![alt text](https://github.com/abuultjens/accessory-sites/blob/master/aln-new.png)






![alt text](https://github.com/abuultjens/accessory-sites/blob/master/snp-sites.png)



![alt text](https://github.com/abuultjens/accessory-sites/blob/master/accessory.png)






