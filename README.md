# Tip and tricks for VCF files

- [Usefull tools](https://github.com/IARC-bioinfo/VCF-tricks#usefull-tools)
- [Use R VariantAnnotation bioconductor package](https://github.com/IARC-bioinfo/VCF-tricks#use-r-variantannotation-bioconductor-package)
- [Manually processing VCF in R](https://github.com/IARC-bioinfo/VCF-tricks#manually-processing-vcf-in-r)
- [Manually processing VCF in bash](https://github.com/IARC-bioinfo/VCF-tricks#manually-processing-vcf-in-bash)

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
- R bioconductor package [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html)
- vcflib [github page](https://github.com/ekg/vcflib)
- vt [wiki](http://genome.sph.umich.edu/wiki/Vt) and [github page](https://github.com/atks/vt)
- bedtools [documentation](http://bedtools.readthedocs.org) and [github page](https://github.com/arq5x/bedtools2)
- PyVCF [github page](https://github.com/jamescasbon/PyVCF)
- VCFtools [webpage](https://vcftools.github.io/) and [github page](https://github.com/vcftools/vcftools). It has been mostly replaced with bcftools but some commands are still only available in VCFtools (in particular [vcf-annotate](https://vcftools.github.io/perl_module.html#vcf-annotate))

## Use R VariantAnnotation bioconductor package

All the commands below assume the package `VariantAnnotation` has been loaded into R using `library(VariantAnnotation)`.

### Replace INFO/DP field with GENO/DP field
```R
vcf <- readVcf("test.vcf", "hg19")
info(vcf)$DP=geno(vcf)$DP
writeVcf(vcf,"test.vcf")
```

### Create a new INFO field

Here it's called `DP_T`and filled with `.` (dot represent missing values in VCF files) but it could be anything you like.
```R
vcf <- readVcf("test.vcf", "hg19")
newInfo <- DataFrame(Number=1, Type="Integer",Description="DP in normal",row.names="DP_N")
info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
info(vcf)$DP_N="."
writeVcf(vcf,"test.vcf")
```

### Split a multi-sample VCF into n-single sample VCFs
```R
vcf_file = "test.vcf"
for (cur_sample in samples(scanVcfHeader(vcf_file))) {
  writeVcf(readVcf(vcf_file, "hg19",ScanVcfParam(sample = cur_sample)),paste(cur_sample,".vcf",sep = ""))
}
```

## Manually processing VCF in R

Look at these functions too: https://github.com/sahilseth/vcfparser

### Extract expected Q-value from a needlestack calling
Once a needlestack calling has been launch, user can want to compute what would be the Q-value of a particular sample at a particular position.
This could help controlling the false negative rate.  
Given a VCF chunk from [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html) bioconductor package, a sample ID and a particular mutation, [this](https://github.com/IARCbioinfo/VCF-tricks/blob/master/code/get_Qvalue.r) function extract the corresponding Q-value.

### Loading a VCF file as a data frame
On Unix systems (Mac or Linux), automatically pipe it with `grep` and `sed` to remove the header.
```R
my_vcf=read.table(pipe("grep -v '^##' test.vcf | sed s/^#//"),stringsAsFactors=F,header=T,sep="\t")
```
On Windows, you can manually remove the header lines (starting with `##`) and also the `#` character from the line containing the column names. After that you can read it using:
```R
my_vcf=read.table("test_noheader.vcf",stringsAsFactors=F,header=T,sep="\t")
```

### Two R functions to extract values from INFO or GENOTYPE fields

Gist: https://gist.github.com/mfoll/a4dfbb92068dc559f130
```R
get_info=function(info,field,num=T) {
  get_single_info=function(single_info,field) {
    grep_res=grep(paste("^",field,"=",sep=""),unlist(strsplit(single_info,";")),value=T)
    if (length(grep_res)>0) strsplit(grep_res,"=")[[1]][2] else NA
  }
  res=unlist(lapply(info,get_single_info,field))
  if (num) as.numeric(res) else res
}

get_genotype=function(genotype,format,field,num=T) {
  get_single_genotype=function(single_genotype,format,field) {
    single_res=unlist(strsplit(single_genotype,":"))[which(unlist(strsplit(format,":"))==field)]
    if (length(single_res)>0) single_res else NA
  }
  res=unlist(lapply(genotype,get_single_genotype,format,field))
  if (num) as.numeric(res) else res
}
```

Both function are vectorized (i.e. you can give them a vector of `INFO` fields or a vector of `GENOTYPE` fields. The genotype field requires that you give the format of the field (for example `"GT:AO:DP"`). In both functions the `field` argument indicates which field you want to extract. By default the result is converted to a numeric value, unless you specify `num=FALSE` when you call the functions.

### Get genotype columns and sample names:
```R
# list of columns containing sample specific data
GT_cols=(which(names(my_vcf)=="FORMAT")+1):ncol(my_vcf)
# extract sample names
SM=names(my_vcf)[GT_cols]
```

### Get number of variants for each position in a VCF file
Following function returns, for each position in the input VCF file, the number of samples having a `QVALUE` higher than the input threshold (`QVAL_thr`, default=50), __i.e__. the number of variants.  

```R
# after reading a VCF with read.table(), see last example
get_number_of_variants <- function(vcf, QVAL_thr = 50){
  unlist(lapply(1:nrow(vcf), function(i) {
    all_QVAL = unlist(lapply((which(colnames(vcf)=="FORMAT")+1) : ncol(vcf),
                             function(id) get_genotype(vcf[i,id], vcf[i,"FORMAT"], field = "QVAL") ))
    sum(all_QVAL>=QVAL_thr)
  }))
}
```

### Using all the above
This assumes that you have a `my_vcf` data frame loaded, the two functions above and the objects `GT_cols` and `SM`.

- Extract the variant type (`TYPE`) from all lines from the INFO field:

  ```R
  get_info(my_vcf$INFO,"TYPE",num=F)
  ```
- Use it to fikter only variants with `TYPE=snv`:

  ```R
  my_vcf[which(get_info(my_vcf$INFO,"TYPE",num=F)=="snv"),]
  ```
- Extract coverage (`DP`) of each sample at a given line (1 here):

  ```R
  get_genotype(my_vcf[1,GT_cols],my_vcf$FORMAT[1],"DP")
  ```

- Extract coverage of all lines of a given sample (`MY_SAMPLE` here):

  ```R
  get_genotype(my_vcf[,"MY_SAMPLE"],my_vcf$FORMAT[1],"DP")
  ```
You can replace `"MY_SAMPLE"` with `SM[1]` to take the first sample without typing manually its name (usefull if you have only one for example).

- Plot the distribution of allelic fraction from the first sample (assuming `AO` and `DP` fields are available in the genotype column):

  ```R
  AO=get_genotype(my_vcf[,SM[1]],my_vcf$FORMAT[1],"AO")
  DP=get_genotype(my_vcf[,SM[1]],my_vcf$FORMAT[1],"DP")
  hist(AO/DP)
  ```

### Built a TSV file, from VCF to extract a particular field from FORMAT

[extract_FORMAT_vcf.r](https://github.com/IARCbioinfo/VCF-tricks/blob/master/code/extract_FORMAT_vcf.r) is a script which extract, for each variant in the VCF, and for each sample, the value of a `field` in `FORMAT`.  
Example of command line:
```
Rscript extract_FORMAT_vcf.r --input_vcf=testepic.vcf.gz --field=DS
```


## Manually processing VCF in bash

### Split a n-samples VCF

This bash script splits a big VCF from n samples into n VCF with file name = sample name (save these lines into big_VCF_to_samples.sh)

Gist : https://gist.github.com/tdelhomme/cb28dec176b55c43e887
```bash
#!/bin/bash

if [ $# -eq 0 ]; then #if no provided parameters
  echo 'usage : ./big_VCF_to_samples.sh <input big VCF> <result folder>'
else
  mkdir -p $2
  IFS= read -a array <<< $(grep "#CHROM" $1 | head -1 | awk '{for(i=10;i<=NF;++i)print $i}')
  samples=${array[0]}
  for i in `seq 1 $(echo "$samples" | wc -w)`;
  do
    bcftools view -s $(echo "$samples" | cut -d" " -f$i) $1 > $2"/"$(echo "$samples" | cut -d" " -f$i).vcf
  done
fi
```

### Generate a sorted and merged BED file from positions in a VCF
```
awk '{ if (!/^#/) print $1"	"$2"	"$2+1}' input.vcf | sort -k1,1 -k2,2n | bedtools merge -i stdin
```

### Compute number of positions in a bed file
```
cat mybedfile.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2+1 }END{print SUM}'
```

### Get the samples ID from a VCF
```
grep "^#CHROM" input.vcf | tr '\t' '\n' | grep -v -E '#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT'
```

### Downsample a VCF
Here a zipped one (assuming the header is less than 10000 lines), from which we randomly extract 1 million lines:
```
(zcat human_9606_b150_GRCh38p7.vcf.gz | head -n 10000 | zgrep ^# ; zgrep -v ^# human_9606_b150_GRCh38p7.vcf.gz | shuf -n 1000000 | LC_ALL=C sort -k1,1V -k2,2n) | gzip > human_9606_b150_GRCh38p7_small.vcf.gz
```
