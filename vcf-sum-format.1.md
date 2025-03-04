% vcf-sum-format(1) Version 0.9 | General Commands

# NAME

**vcf-sum-format** - Sum FORMAT fields of multiple VCF files

# SYNOPSIS

| **vcf-merge-info** \[_options_\] _input.vcf_ \[...\]

# DESCRIPTION

Sums the data in the FORMAT fields from multiple VCF files. Note that this will not permute the sample names, use **vcf-add-sample** to do that, but will blindly add the corresponding samples.

## Options

-h 
: Use VCF headers

-S _sample_,\[...\]
: Rename samples
