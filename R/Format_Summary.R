################################################################## Required Packages

################################################################## Function Excluding Overlap Samples

#' Format Summary
#'
#' Reformat GWAS summary statistics file before overlap_adjustment
#' @importFrom data.table fread
#' @import dplyr
#' @param input_file Path of Summary statistics to be formatted
#' @param Genome_Build Genome build. hg37 or hg38
#' @param Pop Population group (EUR or EAS)
#' @param cols Column names of Summary statistics
#' @param Spcf_n Specify sample size if not stated in column
#' @param Spcf_n_eff Specify effective sample size if not stated in column
#' @param Spcf_n_case Size of overlapped samples if not stated in column
#' @param Spcf_n_ctrl Effective Sample size of summary to be adjusted
#' @param phenotype Effective Sample Size of the overlapped
#' @param allele_flip Sample Size of the all cases
#' @param filter_by_hapmap3 Sample Size of the all controls
#' @param fill_missing Imputation missing columns with given Summary statistics
#' @param hapmap3_only Scope variants into hapmap3 markers only
#' @param minimum.P Minimum pvalue supported : default=9.88e-324
#' @param save_path Path for saving formatted summary statistics (optional)
#' @return Formatted Summary Statistics
#' @examples
#' cs = Check_Sums(system.file('extdata/example/consortium.ss',package = 'OAPRS'),
#' Genome_Build = "hg37", Pop = "eas",
#' cols = c(BETA="beta",Pval="pval",CHR="chrom",POS="pos",REF="ref",ALT="alt",SNP="rsids"),
#' Spcf_n=249625,Spcf_n_case = 50466, Spcf_n_ctrl = 199159)
#' @export
Check_Sums = function(
    input_file,
    Genome_Build,
    Pop,
    cols= c(BETA="BETA",SE="SE",Z=NULL,MAF="MAF",Pval="P",
            CHR="CHR",POS="POS",REF="Allele1",ALT="Allele2",SNP="SNPID",
            n="N",n_eff=NULL,n_case=NULL, n_ctrl=NULL),
    Spcf_n = NULL, Spcf_n_eff = NULL, Spcf_n_case = NULL, Spcf_n_ctrl = NULL,
    phenotype = "binary",
    allele_flip = TRUE,
    filter_by_hapmap3 = TRUE,
    fill_missing = TRUE,
    hapmap3_only=TRUE,
    minimum.P=9.88e-324,
    save_path = NULL
)
{
  BETA<-CHR<-BP<-SNP<-MAF<-POS<-REF<-ALT<-ss_new<-NULL
  n_all<- n_case <- n_ctrl<- n_ov <- NULL
  headers = colnames(fread(input_file,nrows=1))
  #print(headers)
  if (!is.numeric(cols)){
    msg1="Wrongly Specified Column Names: "
    msg2=c()
    for ( i in cols){
      if (!i %in% headers){
        msg2=c(msg2,i)
      }
    }
    if (length(msg2)!=0){stop(paste0(msg1,paste(msg2,collapse=", ")))}
  }
  if ((is.na(cols['BETA'])|is.na(cols['SE'])) & is.na(cols['Z']) & is.na(cols['Pval'])) {
    stop("Please Provide either Pvalue, Z-score or effect size column with SE column")
  }
  if ((is.na(cols['n']))) {
    if(is.null(Spcf_n)){
      stop("Please Provide Sample Size")
    }
  }
  if (phenotype == "binary"){
    if ((is.na(cols['n_case'])|is.na(cols['n_ctrl'])) & is.na(cols['n_eff'])) {
      if(is.null(Spcf_n_case)|is.null(Spcf_n_ctrl) & is.null(Spcf_n_eff)){
        stop("Please Provide either case/control size or effective sample size")
      }
    }
  }else if(phenotype=="quantitative"){
    Spcf_n_eff=Spcf_n; cols['n_eff']=cols['n']
  }else{
    stop("Please Specify Phenotype: 'binary' or 'quantitative'")
  }
  if (!is.na(Genome_Build) & Genome_Build%in%c("hg37","hg38")){
    msg = paste0("Checking Summary Statistics in ",Genome_Build)
  }else{
    stop("Please Specify Genome_Build: 'hg37' or 'hg38'")
  }

  if (!Pop %in% c("eur","eas")){
    message("OAPRS supports LD independent blocks for only EUR and EAS Population")
    stop("Please Specify Population Group: 'eur' or 'eas'")
  }
  if (filter_by_hapmap3){
    hm = fread(system.file(paste0('extdata/SNPInfo/snpinfo_1kg_hm3_',Pop,".gz"),package = "OAPRS"))
  }
  message(paste0(msg," in ",toupper(Pop), " Population"))
  if (is.numeric(cols)){tmp=headers[cols]; names(tmp) = tmp; cols = tmp}
  print(cols)
  ss = fread(input_file,fill=TRUE) %>% select(all_of(cols))
  if (is.na(cols['n'])){ss$n = Spcf_n}
  if(phenotype=="continuous"){ss$n_eff = ss$n}
  if(phenotype=="binary"){
  if (is.na(cols['n_case'])){ss$n_case = Spcf_n_case}
  if (is.na(cols['n_ctrl'])){ss$n_ctrl = Spcf_n_ctrl}
  if (!is.na(cols['n_eff'])){
    n_eff = ss$n_eff
  }else if(!is.null(Spcf_n_eff)){ss$n_eff = Spcf_n_eff
  }else{
    ss=ss%>%mutate(n_eff = 4/(1/n_case+1/n_ctrl))
  }
  }
  if (is.na(cols['MAF']) & is.na(cols['SNP'])){ss = left_join(ss, hm %>% select(CHR,BP,SNP,MAF), by = c(setNames('BP', 'POS'),setNames('CHR', 'CHR')))
  }else if(is.na(cols['MAF'])){ss = left_join(ss, hm %>% select(CHR,BP,MAF), by = c(setNames('BP', 'POS'),setNames('CHR', 'CHR')))}
  if (is.na(cols['Z'])&is.na(cols['Pval'])){tmp = ss$SE; tmp2 = ss$BETA; ss$Z=tmp2/tmp}
  if (is.na(cols['Z'])&!is.na(cols['Pval'])){tmp = ss$Pval; tmp2 = ss$BETA; ss$Z = PtoZ(tmp,tmp2)}
  if (is.na(cols['Pval'])){Z = ss$Z;ss = ss %>% mutate(Pval = 2*pnorm(Z,lower.tail = F))}
  if (is.na(cols['SE'])&!is.na(cols['BETA'])){maf = ss$MAF; ss=ss%>%mutate(SE = BETA/Z) }
  if (is.na(cols['SE'])&is.na(cols['BETA'])){maf = ss$MAF; ss=ss%>%mutate(SE = 1/sqrt(2*MAF*(1-MAF)*(n_eff + Z^2))) }
  if (is.na(cols['BETA'])){Z = ss$Z;SE =ss$SE; ss=ss%>%mutate(Beta = Z*SE)}
  ss= ss%>%mutate(CHR_POS = paste0(CHR,":",POS,"_",REF,"/",ALT))
  if (all(ss$Pval<=0)) {message("Assuming Pvalue as log transformed: "); ss$P = exp(ss$P)}
  print(head(ss))
  maf = ss$MAF; SE = ss$SE; Pval = ss$Pval ; Z = ss$Z; n_eff=ss$n_eff
  err_se = SE<=0
  miss_se = is.na(SE)
  message(paste0("Number of Markers with SE errors: ",sum(err_se),"  Number of Markers with SE missing: ",sum(miss_se)))
  err_P = Pval<=0
  miss_P = is.na(Pval)
  message(paste0("Number of Markers with Pval errors: ",sum(err_P),"  Number of Markers with Pval missing: ",sum(miss_P)))
  if (fill_missing){message("Imputing missing..")
    SE_new=1/sqrt(2*maf*(1-maf)*((n_eff) + Z^2)); ss$SE[miss_se|err_se] = ss_new[miss_se|err_se]
    Pval_new = 2*pnorm(Z,lower.tail = F); ss$Pval[err_P|miss_P] = Pval_new[err_P|miss_P]
  }
  err_maf = maf<0 | maf>1
  miss_maf = is.na(maf)
  message(paste0("Number of Markers with MAF errors: ",sum(err_maf),"  Number of Markers with MAF missing: ",sum(miss_maf)))
  ss = ss %>% filter(!(err_maf|miss_maf))
  if(!is.null(save_path)){write.table(ss,save_path,col.names = T, row.names = F, quote = F)}
  return(ss)
}
