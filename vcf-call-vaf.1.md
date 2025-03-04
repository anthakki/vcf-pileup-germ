% vcf-call-vaf(1) Version 0.9 | General Commands

# NAME

**vcf-call-vaf** - Flag high/low relative variant allele frequency

# SYNOPSIS

| **vcf-call-vaf** \[_options_\] _input.vcf_

# DESCRIPTION

Flags variants with high/low variant allele frequency (**HIVAF** and **LOVAF** in the VCF INFO, respectively) in the tumor relative to the normal. The VCF needs to be first piled up with **vcf-pileup**.

## Options

-h 
: Use VCF headers

-v _vaf_or_
: Threshold for tumor-to-normal VAF

-s _sample_
: Mark sample _sample_ as normal (germline only)

-S _sample_
: Mark sample _sample_ as tumor
