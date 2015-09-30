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
- R bioconductor package [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html)
- vcflib [github page](https://github.com/ekg/vcflib)
- vt [wiki](http://genome.sph.umich.edu/wiki/Vt) and [github page](https://github.com/atks/vt)
- bedtools [documentation](http://bedtools.readthedocs.org) and [github page](https://github.com/arq5x/bedtools2)

## R Tip and tricks

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
GT_cols=(which(names(data)=="FORMAT")+1):ncol(data)
# extract sample names
SM=names(data[,GT_cols])
```

### Using all the above
This assumes that you have a `my_vcf` data frame loaded and the two function above.

- Extract the variant type (`TYPE`) from all lines from the INFO field:
  
  ```R
  get_info(my_vcf$INFO,"TYPE",num=F)
  ```
- Use it to fikter only variants with `TYPE=snv`:

  ```R
  data[which(get_info(my_vcf$INFO,"TYPE",num=F)=="snv"),]
  ```
- Extract coverage (`DP`) of each sample at a given line (1 here):

  ```R
  get_genotype(my_vcf[1,GT_cols],my_vcf$FORMAT[1],"DP")
  ```
  
- Extract coverage of all lines of a given sample (`MY_SAMPLE` here):

  ```R
  get_genotype(my_vcf[,"MY_SAMPLE"],my_vcf$FORMAT[1],"DP")
  ```
You can replace `"MY_SAMPLE"` with `SM[1]` to take the first sample without typing manually its name.
