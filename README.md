
# Tools for post-hoc VCF germline pileup/calling 

## Prerequisites

To build the C++ programs, a C++ 11 and a C 99 compiler, GNU `make`, `pkg-config`, and [`htslib`](https://github.com/samtools/htslib) (of SAMtools fame) and [`zlib`](https://zlib.net/) are needed. To run some of the tools, Python 3 is needed.

To build the binaries, simply run:
```
	make
```

The software has been tested on macOS 14 and Amazon Linux 2.

## Usage

The [`vcf-pileup`](vcf-pileup.1.md) program allows the pileup of read counts of specific variants, as specified by a VCF file, from the original BAMs. For example, to pileup existing somatic calls from `somatic.vcf` into `pileup.vcf` from the BAMs `tumor.bam` and `normal.bam`:
```
	./vcf-pileup -h somatic.vcf tumor.bam normal.bam >pileup.vcf
```
The [`vcf-call-germ`](vcf-call-germ.1.md) program allows genotyping based on the pileup, and can flag germline variants in the resulting VCF. This runs the analysis with default parameters, adding information to `flagged.vcf`:
```
	./vcf-call-germ -h -S tumor pileup.vcf >flagged.vcf
```

This is likely all that is needed. See the manual pages above for more options and details. See [`gtbeta.c`](gtbeta.c) for the model.

The [`vcf-call-fs`](vcf-call-fs.1.md) program can flag strand bias based on the pileup using standard metrics [`FS`](https://gatk.zendesk.com/hc/en-us/articles/360036896952-FisherStrand) and [`SOR`](https://gatk.zendesk.com/hc/en-us/articles/360037267551-AS-StrandOddsRatio) instead of Bayesian calling:
```
	./vcf-call-fs -h pileup.vcf >flagged.vcf
```

The [`vcf-call-vaf`](vcf-call-vaf.1.md) is a simple script to flag variants by tumor-to-normal VAF based on the pileup:
```
	./vcf-call-vaf -h -S tumor pileup.vcf >flagged.vcf
```

The scripts [`vcf-add-samples`](vcf-add-samples.1.md), [`vcf-merge-info`](vcf-merge-info.1.md), [`vcf-sum-format`](vcf-sum-format.1.md), and 
[`vcf-header`](vcf-header.1.md) might be useful in manipulating the resulting pileup VCFs. [`bcftools`](https://github.com/samtools/bcftools) [`filter`](https://samtools.github.io/bcftools/bcftools.html#filter) is useful in filtering the flagged variants.

## Copying

The files are distributed under the [MIT license](LICENSE.txt). Copyright (c) 2023 Antti Hakkinen, except portions of [`jacobi_rule.c`](jacobi_rule.c).
