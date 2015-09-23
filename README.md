# Tip and tricks for VCF files

## Usefull tools
### Samtools organisation and repositories
- File format [specification](http://samtools.github.io/hts-specs/)
- Bcftools [github page](https://github.com/samtools/bcftools)
- Bcftools [webpage](http://samtools.github.io/bcftools/)
- Bcftools [man page](http://samtools.github.io/bcftools/bcftools.html)

Compilation (from [here](http://samtools.github.io/bcftools/)):
```bash
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
git clone --branch=develop git://github.com/samtools/samtools.git
cd bcftools; make
cd ../samtools; make
cd ../htslib; make
```

### Other tools
- R bioconductor package [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.ht
- vcflib [github page](https://github.com/ekg/vcflib)
- vt [wiki](http://genome.sph.umich.edu/wiki/Vt) and [github page](https://github.com/atks/vt)
- bedtools [documentation](http://bedtools.readthedocs.org) and [github page](https://github.com/arq5x/bedtools2)

## Tip and tricks

### Loading a VCF file as a data frame

```R
my_vcf=read.table(pipe("grep -v '^##' test.vcf | sed s/^#//"),stringsAsFactors=F,header=T,sep="\t")
```

Note that the header is ignored.

### Two R functions to extract values from INFO or GENOTYPE fields 

<script src="https://gist.github.com/mfoll/a4dfbb92068dc559f130.js"></script>

Usage example:
