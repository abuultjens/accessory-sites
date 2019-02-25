# accessory-sites
A script to extract invariant accessory sites from whole genome alignments

# Author
Andrew Buultjens

# Synopsis
Accessory sites are positions in an alignment in which one or more taxa have missing information that is generally denoted by 'N'. Such missing information might be due to a number of interesting reasons such as a deletion, coverage below a specific snp calling threshold or mixed alleles. Accessory sites not containing snps, hereon referred to as 'invariant accessory sites' are often overlooked in comparative genomic studies, however they may add an additional level of information alongside that of snps in core and accessory sites.

# Usage
$ sh accessory-sites.sh [INFILE.fa] [PREFIX] [OUTFILE_FORMAT] [OUTFILE_DATA]  

* note that the command must have the arguments in the exact order as specified above

# Options
* [OUTFILE_FORMAT] can be either 'table' or 'fasta'   
* [OUTFILE_DATA] can be either 'all_sites' or 'only_invariant_accessory_sites'

# How it works
accessory-sites tricks the incredibly fast snp-sites into thinking that invariant accessory positions actually contain snps so that they are extracted from alignments. This is useful as the fake snp sites are then converted back to being invariant accessory sites and concatenated to the original variant core and accessory sites extracted by snp-sites.

Basically, snp-sites is first run with the original alignment and the sites extracted are recorded. All 'N' in the alignment are then replaced with 'A', effectively introducing a fake snp, and this modified alignment is run through snp-sites with any new sites compared to the original run identified. Novel sites are then extracted, 'A' replaced with the original 'N' and the invariant accessory positions are added to the outfile. This process is repeated for 'G', 'C' and 'T' so that all invariant accessory sites are obtained.



# Dependencies
* vcf-tools     
* bedtools     
* snp-sites    

# Example

![alt text](https://github.com/abuultjens/accessory-sites/blob/master/aln-new.png)






![alt text](https://github.com/abuultjens/accessory-sites/blob/master/snp-sites_output.png)



![alt text](https://github.com/abuultjens/accessory-sites/blob/master/output.png)






