library(VariantAnnotation)

args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$input_vcf))            {stop("no input VCF file")} else {input_vcf = args$input_vcf}
if(is.null(args$field))                {stop("no FORMAT field to query")} else {field=args$field}
if(is.null(args$out_txt))              {out_txt = paste(gsub(".vcf.gz","",input_vcf),"_",field,"_extract.txt",sep="")} else {out_txt=args$out_txt}
if(is.null(args$nb_mut))               {nb_mut = 1000000} else {nb_mut = as.numeric(args$nb_mut)}

vcf <- open(VcfFile(input_vcf,  yieldSize=nb_mut))
vcf_chunk = readVcf(vcf, "hg19")
`#CHR` = as.character(seqnames(rowRanges(vcf_chunk)))
`#START` = start(ranges(rowRanges(vcf_chunk,"seqnames")))
`#REF` = as.character(rowRanges(vcf_chunk)$REF)
`#ALT` = as.character(unlist(rowRanges(vcf_chunk)$ALT))

res_matrix = rbind(`#CHR`, `#START`, `#REF`, `#ALT`, t(geno(vcf_chunk, field)))

write.table(data.frame("." = rownames(res_matrix), res_matrix, check.names = F), file = out_txt, row.names=FALSE, sep="\t", quote = F)
