
#' P-value to z-score
#'
#' P-value to z-score transformation
#' @param P P-values to convert
#' @param bet Signed vector of coefficients
#' @export
PtoZ = function(P,bet){
  return(abs(qnorm(P/2))*sign(bet))
}

#' Overlap Adjustment
#'
#' Overlap Adjustment from a GWAS summary with summary of Overlapped samples
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @import dplyr
#' @param file_all object/Path of Summary statistics to be adjusted
#' @param file_ov object/Path of Summary statistics of Overlapped samples
#' @param save_path Output path of adjusted summary statistics
#' @param phenotype Phenotype if binary : "binary", continuous : "continuous"
#' @param dropna Drop markers that has standard error NA's?
#' @return Adjusted summary
#' @examples
#' print('adj_ss = exclude_overlap(cs,ts,"adj.txt",phenotype="binary")')
#' @export
exclude_overlap = function
(
  file_all,
  file_ov,
  save_path=NULL,
  phenotype,
  dropna=T
)
{
  AF_Allele2_all<-BETA_IVW<-BETA_all<- BETA_ov<- P_Part<- SE_IVW <- NULL
  SE_Part<- SE_all<- SE_ov<- case_af <- chr <- control_af <- NULL
  MAF_all <- n_all<- n_case <- n_ctrl<- n_ov <- NULL
  # File read
  if(is.character(file_all)){file_all = fread(file_all,fill=TRUE)}
  if(is.character(file_ov)){file_ov = fread(file_ov,fill=TRUE)}
  ss_all = file_all
  ss_ov = file_ov
  check_cols = c("BETA",  "Pval","CHR","POS","REF","ALT","SNP","n","n_case", "n_ctrl","n_eff","MAF","Z","SE","CHR_POS")
  if(!all(check_cols%in%colnames(ss_all))){stop('Please Check Format of Summary')}

  # Subset only intersection of chromosal position of each GWAS
  ss_ov$REF = toupper(ss_ov$REF);ss_ov$ALT = toupper(ss_ov$ALT);ss_all$REF = toupper(ss_all$REF);ss_all$ALT = toupper(ss_all$ALT)
  ss_tmp = inner_join(ss_ov,ss_all, by = c("CHR","POS","REF","ALT"),suffix = c("_ov", "_all"))
  ss_tmp2 = inner_join(ss_ov,ss_all, by = c("CHR","POS","REF"="ALT","ALT"="REF"),suffix = c("_ov", "_all"))
  ss_tmp2 = ss_tmp2%>%mutate(BETA_all=(-1)*BETA_all,MAF_all=1-MAF_all)
  ss_tmp = rbind(ss_tmp,ss_tmp2)
  ss_tmp$Pval_all = as.numeric(ss_tmp$Pval_all)
  ss_tmp$Pval_ov = as.numeric(ss_tmp$Pval_ov)
  print("Number of variant of all, overlapped, subtracted samples")
  print(c(nrow(ss_all),nrow(ss_ov),nrow(ss_tmp)))

  # Sample Size Calculation
  est_all = ss_tmp$BETA_all; se_all = ss_tmp$SE_all ;
  est_ov = ss_tmp$BETA_ov; se_ov =ss_tmp$SE_ov ;
  n_all=ss_tmp$n_all; n_ov=ss_tmp$n_ov
  n_exc = ss_tmp$n_all-ss_tmp$n_ov

  # Get Minor Allele Frequency
  maf_all = ss_tmp$MAF_all
  maf_ov = ss_tmp$MAF_ov
  maf=(maf_all*n_all - maf_ov*n_ov)/n_exc

  # Inverse Variance Weighting
  ss_tmp = ss_tmp%>%mutate(SE_IVW = sqrt(1/((1/SE_all^2)-(1/SE_ov^2))))
  ss_tmp = ss_tmp%>%mutate(BETA_IVW = (BETA_all/SE_all^2-BETA_ov/SE_ov^2)*SE_IVW^2)
  ss_tmp = ss_tmp%>%mutate(Pval_IVW = exp(pchisq((BETA_IVW^2)/SE_IVW^2,df=1,lower.tail=F,log.p=T)))

  # Reversed z-score
    eff_all =ss_tmp$n_eff_all; eff_ov = ss_tmp$n_eff_ov
    w_all = sqrt(eff_all); w_ov = sqrt(eff_ov); w_RZ = sqrt(eff_all-eff_ov)
    z_all = PtoZ(ss_tmp$Pval_all,ss_tmp$BETA_all); z_ov = PtoZ(ss_tmp$Pval_ov,ss_tmp$BETA_ov)
    z_RZ = z_all*w_all - z_ov*w_ov ; z_RZ = z_RZ/w_RZ
    SE_RZ = 1/sqrt(2*maf_all*(1-maf_all)*((eff_all-eff_ov) + z_RZ^2))
    BETA_RZ =  z_RZ * SE_RZ
    sign_RZ = sign(z_RZ); Pval_RZ = exp(pnorm(abs(z_RZ),lower.tail=F,log.p=T))*2
  # Generate Output GWAS Dataframe
  new_ss = data.frame(
        CHR = ss_tmp$CHR, POS = ss_tmp$POS, ALT = ss_tmp$ALT,
        REF = ss_tmp$REF, MAF = maf, MAF_all = maf_all, N = n_exc, N_all = ss_tmp$n_all, N_ov = ss_tmp$n_ov,
        BETA_IVW = ss_tmp$BETA_IVW,SE_IVW = ss_tmp$SE_IVW,P_IVW = ss_tmp$Pval_IVW,
        BETA_all = ss_tmp$BETA_all,P_all = ss_tmp$Pval_all,SE_all = ss_tmp$SE_all,
        BETA_ov = ss_tmp$BETA_ov,P_ov = ss_tmp$Pval_ov, SE_ov = ss_tmp$SE_ov, SNP = ss_tmp$SNP_ov, CHR_POS=ss_tmp$CHR_POS_all)
  new_ss = new_ss %>% mutate(sign_RZ = sign_RZ, P_RZ = Pval_RZ, z_RZ = z_RZ, z_all = z_all, z_ov = z_ov, BETA_RZ = BETA_RZ, SE_RZ = SE_RZ)
  if (!dropna){
  new_ss$SE_IVW[is.na(new_ss$SE_IVW)] = max(ss_tmp$SE_all[is.na(new_ss$SE_IVW)],ss_tmp$SE_ov[is.na(new_ss$SE_IVW)])
  new_ss$SE_Part[is.na(new_ss$SE_Part)] = max(ss_tmp$SE_all[is.na(new_ss$SE_Part)],ss_tmp$SE_ov[is.na(new_ss$SE_Part)])
  }
  new_ss = new_ss %>% filter(!is.na(SE_IVW))
  # Check header and size of new overlap excluded GWAS
  print("New Summary")
  print(head(new_ss))
  print("Dimension of new summary")
  print(dim(new_ss))
  if(!is.null(save_path)){write.table(new_ss,save_path,col.names = T, row.names = F, quote = F)}
  return(new_ss)
}

