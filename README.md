# accessory-sites
A script to extract invariant accessory sites (accessory genome positions with no snps) from whole genome alignments

# Author
Andrew Buultjens

# Synopsis
Accessory sites are positions in an alignment in which one or more taxa have missing information that is generally denoted by 'N'. Such missing information might be due to a number of biologically interesting reasons such as a deletion, coverage below a specific snp calling threshold or mixed alleles. Accessory sites not containing snps, hereon referred to as 'invariant accessory sites' are often overlooked in comparative genomic studies, however they may add an additional level of information alongside that of snps in core and accessory sites. snp-sites (https://github.com/sanger-pathogens/snp-sites) is a tool that rapidly extracts core and accessory snp sites from whole genome alignments, however it currently does not allow for the extraction of invariant accessory sites. accessory-sites builds upon the speed of snp-sites to extract invariant accessory sites from whole genome alignments.

# Usage
$ sh accessory-sites.sh [INFILE.fa] [PREFIX] [OUTFILE_DATA]  

* note that the command options must be in the exact order as specified above as they are treated as positional arguments

# Options
* [OUTFILE_DATA] can be either 'all_sites' or 'only_invariant_accessory_sites'

# Help
```
sh accessory-sites.sh help
```

# How it works
accessory-sites tricks the incredibly fast snp-sites into thinking that invariant accessory positions actually contain snps so that they are extracted from alignments. This is useful as the fake snp sites are then converted back to being invariant accessory sites and concatenated to the original variant core and accessory sites extracted by snp-sites.

Basically, snp-sites is first run with the original alignment and the sites extracted are recorded. All 'N' in the alignment are then replaced with 'A', effectively introducing a fake snp, and this modified alignment is run through snp-sites with any new sites compared to the original run identified. Novel sites are then extracted, 'A' replaced with the original 'N' and the invariant accessory positions are added to the outfile. This process is repeated for 'G', 'C' and 'T' so that all invariant accessory sites are obtained.

# Dependencies
* vcf-tools     
* bedtools     
* snp-sites    

# Example

**original alignment:**   
```
cat test.fa
>TAXA_A
ATTN
>TAXA_B
AGTA
>TAXA_C
AGAA
>TAXA_D
AGNN
```
**The alignment contains:**   
Site 1 = core position no snp (no information)   
Site 2 = core position with snp (useful information)   
Site 3 = accessory position with snp (useful information)   
Site 4 = accessory position with no snp (useful information)   

**output from snp-sites:**   
```
snp-sites -v test.fa | vcf-to-tab | tr -d '/' | tr '*' 'N' | cut -f 1,2,4- 
#CHROM	POS	TAXA_A	TAXA_B	TAXA_C	TAXA_D
1	2	T	G	G	G
1	3	T	T	A	N
```
**sites kept:**  
Site 2 = core position with snp  
Site 3 = accessory position with snp  

**sites removed:**  
Site 1 = core position no snp  
Site 4 = invariant accessory position (accessory position with no snp)  

# Example 1: run accessory-sites to extract invariant accessory sites and all snp sites:

```
sh accessory-sites.sh test.fa OUT all_sites
```
inspect fasta output
```
cat OUT.fa
>TAXA_A
TTN
>TAXA_B
GTA
>TAXA_C
GAA
>TAXA_D
GNN
```
inspect table output
```
cat OUT.tab
#CHROM	POS	TAXA_A	TAXA_B	TAXA_C	TAXA_D
1	2	T	G	G	G
1	3	T	T	A	N
1	4	N	A	A	N
```
**sites kept:**  
Site 2 = core position with snp  
Site 3 = accessory position with snp  
Site 4 = invariant accessory position (accessory position with no snp)  

**sites removed:**  
Site 1 = core position no snp  

# Example 2: run accessory-sites to extract only invariant accessory sites:
```
sh accessory-sites.sh test.fa OUT only_invariant_accessory_sites
```
inspect fasta output
```
cat OUT.fa
>TAXA_A
N
>TAXA_B
A
>TAXA_C
A
>TAXA_D
N
```
inspect table output
```
cat OUT.tab
#CHROM	POS	TAXA_A	TAXA_B	TAXA_C	TAXA_D
1	4	N	A	A	N
```
**sites kept:**  
Site 4 = invariant accessory position (accessory position with no snp)  

**sites removed:**   
Site 1 = core position no snp  
Site 2 = core position with snp   
Site 3 = accessory position with snp  

# Removing a reference sequence from a multi-fasta file
To remove a specific entry from a multi-fasta file (eg. 'Reference') use ref-remover:
https://github.com/abuultjens/ref-remover

