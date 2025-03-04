% vcf-add-samples(1) Version 0.9 | General Commands

# NAME

**vcf-add-samples** - Add (or remove) samples to a VCF file

# SYNOPSIS

| **vcf-add-samples** \[_options_\] _input.vcf_ _sample_ \[...\]

# DESCRIPTION

Generates a VCF where the samples are _sample_ \[...\]. This implies that the
existing samples are permuted, the extra samples are removed, and empty samples
are synthesized for the missing. When synthesizing new samples, the script
tries to infer correct field format from the existing data.

Alternatively, the **-v** option can be used to drop samples by name.

## Options

-h 
: Use VCF headers

-v 
: Use all existing samples but those in _sample_ \[...\]
