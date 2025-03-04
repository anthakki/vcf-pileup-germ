% vcf-call-germ(1) Version 0.9 | General Commands

# NAME

**vcf-call-germ** - Call genotypes and flag germline variants

# SYNOPSIS

| **vcf-call-germ** \[_options_\] _input.vcf_

# DESCRIPTION

Calls genotypes using the **gtbeta22** model and flags germline variants. The VCF needs to be first piled up with **vcf-pileup**. A variant is flagged germline (**GERMLINE** in the VCF INFO field) if it so more probably than not.

The samples can be marked normal (germline only) or tumor (mixed somatic and germline) with **-s** and **-S** respectively, with the default being the opposite of what was last used, or tumor otherwise.

The model involves a common error, whose prior can be controlled by **-e**, **-ew**, modeling end-to-end calling error (as opposed to read quality) at the site; sample specific purity, or mixing, rates, with uniform prior; and biallelic genotypes (0/0, 0/1, 1/1), with uniform prior. The genotypes (**GT**) and genotype log-likelihoods (**GL**) are reported in the VCF FORMAT data for the corresponding sample. Full Bayesian solution is exactly integrated (unless **-at** is used), no Monte Carlo is used, which means the caller can be both fast and accurate.

If **-b** is used, the forward and backward strands are modeled separately as well as together and potential strand bias in the variant is flagged as well (**STRANDBIAS** in the VCF INFO field).

## Options

-h 
: Use VCF headers

-e _em_
: Set error prior mean (default: 0.1)

-ew _ew_
: Set error prior weight (default: 1/_em_)

-b
: Call with bistrand models and flag strand bias as well (default: false)

-Z
: Condition the variant to be somatic (default: false). Can be useful if the variants were already filtered to be somatic with high confidence, but the germline status is questionable.

-s _sample_
: Mark sample _sample_ as normal (germline only)

-S _sample_
: Mark sample _sample_ as tumor

-at _at_
: Approximation taper for numerical integration (default: exact). Models of order _at_ or less will be integrated exactly, but higher orders are logarithmically tapered. Typically, this is set to a multiple of your sequencing depth (e.g. 500 for 100x average depth) such that regular sites are computed exactly, but pathological sites (e.g. >>10,000 reads) can be approximately computed roughly at the same rate.
