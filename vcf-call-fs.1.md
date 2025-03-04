% vcf-call-fs(1) Version 0.9 | General Commands

# NAME

**vcf-call-fs** - Compute Fisher Strand and Strand Odds Ratio scores

# SYNOPSIS

| **vcf-call-fs** \[_options_\] _input.vcf_

# DESCRIPTION

Computes the VCF INFO fields with Fisher Strand (**FS**) and Strand Odds Ratio (**SOR**) scores. The VCF needs to be first piled up with **vcf-pileup**.

The FS score is the p-value for a Fisher's exact test, in negative decibels, with the null hypothesis that the strand specific reads are independent from the reads for and against the mutation.

The SOR score is a logarithmic scaled symmetrical odds ratio estimate of dubious statistical grounding computed by GATK.

## Options

-h 
: Use VCF headers
