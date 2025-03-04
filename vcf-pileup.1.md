% vcf-pileup(1) Version 0.9 | General Commands

# NAME

**vcf-pileup** - Pileup specific variants back from the BAMs

# SYNOPSIS

| **vcf-pileup** \[_options_\] _input.vcf_ _input.bam_ \[...\]

# DESCRIPTION

Goes through a list of variants specified by the VCF and piles them up from the provided BAM files.

The extracted attributes: read depth uniquely matching each allele (**AD**), extra number of reads ambiguously matching alleles (**XD**), total read depth passing quality (**DP**), and the corresponding statistics for the forward strand only (**ADF**, **XDF**, and **DPF**, respectively) are added to the VCF.

## Options

-h 
: Use VCF headers

-6
: Assume the quality is in Illumina1.3+ encoding 

-q _mq_
: Minimum mapping quality for an alignment (default: 0)

-Q _bq_
: Minimum base quality for a base (default: 13) 

-C _cache.vcf_
: Use _cache.vcf_ as a read/write cache. This can be useful in a pipeline.
