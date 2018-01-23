# should build a vcf_chunk in a first step, as the following:
# vcf_file=open(VcfFile(vcf,  yieldSize=10000))
# vcf_chunk = readVcf(vcf_file, "hg19")

# mut is in the form: chr:start_ref/alt
# sm corresponds to one sample as it is defined in the VCF file

toQvalue <- function(vcf_chunk, sm, mut, VAF_coeff=1){
  id_sm = which(colnames(geno(vcf_chunk,"GT"))==sm)
  id_mut = which(rownames(geno(vcf_chunk,"GT"))==mut)
  if(isEmpty(id_mut)) return(0)
  mu_est = info(vcf_chunk[id_mut,])$ERR
  sig_est = info(vcf_chunk[id_mut,])$SIG
  all_DP = unlist(as.list(geno(vcf_chunk[id_mut,],"DP")))
  all_AO = unlist(as.list(geno(vcf_chunk[id_mut,],"AO")))
  DPi = all_DP[id_sm]
  AOi = all_AO[id_sm] * VAF_coeff

  if(AOi>DPi) return(0)
  res = unlist(-10*log10(p.adjust((dnbinom(c(all_AO,AOi),size=1/sig_est,mu=mu_est*c(all_DP,DPi)) +
                               pnbinom(c(all_AO,AOi),size=1/sig_est,mu=mu_est*c(all_DP,DPi),lower.tail = F)))
                   [length(all_DP)+1]))
  if(res>1000){return(1000)} else {return(res)}
}
