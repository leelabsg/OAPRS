################################################################## Required Packages
library(data.table)
library(dplyr)
library(stringr)
library(data.table)
library(parallel)
library(pROC)
library(DescTools)
################################################################## Function Excluding Overlap Samples
# Arguments
# 1. Input files consist of total gwas summary and overlapped sample GWAS summary (file_all, file_ov)
# 2. Specify output file path (output_file)
# 3. Specify column names indicating sample size, effect size, standard error, MAF, Chromosome, Position, Allele1, Allele2 of input, respectively (col1, col2)
#### Or Alternative column names indicating sample size, effect size, standard error, MAF, chromosomal position
# 4. Sample size of each GWAS
## If sample size on each variant is different, indicate the column name
# 5. Either Effective sample sizes / case control size needs to be provided
# 6. Phenotype can be either binary or continuous
# 7. Whether to drop markers that has standard error NA's (otherwise use max standard error between overlapped and all)

# Output
# 0. Summary statistics of All / Overlapped/ IVW/ Partitioned / Z-scoring
# 1. Effect size, Standard Error, Sample size, P-value
# 2. xty from partitioning Method
# 3. z scores and effective sample sizes

# Transformation of Pvalue into Z-score
# chr = sub(".*:^_.*", "", rs_number),pos=gsub(".*:(.*)\\_.*", "\\1", rs_number)

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
#' @param file_all Path of Summary statistics to be adjusted
#' @param file_ov Path of Summary statistics of Overlapped samples
#' @param output_file Output path of adjusted summary statistics
#' @param col_all Column names of Summary statistics to be adjusted
#' @param col_ov Column names of Summary statistics of Overlapped samples
#' @param n_all Sample size of summary to be adjusted
#' @param n_ov Size of overlapped samples
#' @param eff_all Effective Sample size of summary to be adjusted
#' @param eff_ov Effective Sample Size of the overlapped
#' @param n_case_all Sample Size of the all cases
#' @param n_ctrl_all Sample Size of the all controls
#' @param n_case_ov Sample Size of the overlapped cases
#' @param n_ctrl_ov Sample Size of the overlapped controls
#' @param phenotype Phenotype if binary : "binary", continuous : "continuous"
#' @param dropna Drop markers that has standard error NA's?
#' @param MAF external maf information path
#' @param rsid RS snp Marker column
#' @return Adjusted summary
#' @examples
#' exclude_overlap("all_new.txt","/media/leelabsg-storage0/seokho/UKBB/overlapGWAS/Agen_T2D/KoGES_DM.tsv.gz","out_new_1004.txt",n_all = 433540,n_ov = 56918,
#' col_all = c("Beta","SE","EAF","Chr","Pos","Allele1","Allele2","P"),
#' eff_all =  211793, eff_ov = 10879.08,
#' col_ov = c("beta","sebeta","maf","chrom","pos","ref","alt","pval"),phenotype = "binary",dropna=T,
#' MAF = '/media/leelabsg-storage0/seokho/hapmap3/snpinfo_1kg_hm3');
#' temp2 <- F_to_C( c(50, 63, 23) );
#' @export
exclude_overlap = function
(
  file_all,file_ov,output_file, # file paths of Total GWAS summary, overlapped sample GWAS summary, and output
  col_all = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2","P"), # Column names of Total GWAS
  col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2","P"), # Column names of overlapped sample GWAS
  n_all, n_ov, # sample size of Total GWAS and overlapped sample GWAS
  eff_all=NULL, eff_ov=NULL, # effective sample size
  n_case_all = NULL, n_ctrl_all =NULL, n_case_ov = NULL, n_ctrl_ov = NULL, # effective sample size with case/control size using Metal
  phenotype, # Phenotype if binary : "binary", continuous : "continuous"
  dropna=T, # Drop markers that has standard error NA's (otherwise use max standard error between overlapped and all)
  MAF = 'snpinfo_1kg_hm3', # external maf information path,
  rsid = NULL
)
{
  # File read
  ss_all = fread(file_all,fill=TRUE)
  ss_ov = fread(file_ov,fill=TRUE)


  # Handle MAF
  if (!col_all[3]%in%colnames(ss_all)){
  maf_tmp = fread(MAF)
  ss_all = left_join(ss_all, maf_tmp, by = c(setNames("CHR",col_all[4]),setNames("SNP", col_all[5])))
  col_all[3] = "MAF"
  ss_all = ss_all %>% filter(!is.na(MAF))
  }
  if (!col_ov[3]%in%colnames(ss_ov) & !is.null(MAF)){
  maf_tmp = fread(MAF)
  ss_ov = left_join(ss_ov, maf_tmp, by = c(setNames("CHR",col_ov[4]),setNames("BP", col_ov[5])))
  col_ov[3] = "MAF"
  ss_ov = ss_ov %>% filter(!is.na(MAF))
  }
  if (!col_ov[3]%in%colnames(ss_ov)){ss_ov = ss_ov %>% mutate(af = (control_af*67127+case_af*5083)/(67127+5083))}
  # Column selection
  ss_all_old = ss_all
  ss_ov_old = ss_ov
  ss_all = ss_all%>%select(all_of(col_all))
  ss_ov = ss_ov%>%select(all_of(col_ov))

  # Metal
    if(!is.null(n_case_all) & !is.null(n_ctrl_all)){
        eff_all = 4/(1/n_case_all+1/n_ctrl_all)
    }
    if(!is.null(n_case_ov) & !is.null(n_ctrl_ov)){
        eff_ov = 4/(1/n_case_ov+1/n_ctrl_ov)
    }

  # Build chromosomal position and its allele to secure uniqueness
  colnames(ss_all) = colnames(ss_ov)= c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2","P")
  for (i in c("n_all","n_ov","eff_all","eff_ov","n_case_all","n_ctrl_all","n_case_ov","n_ctrl_ov")){
    r=get(i)
    if (is.character(r)&grepl("ov", i, fixed = TRUE)){
    ss_ov[[r]]=ss_ov_old[[r]]
    }
    if (is.character(r)&grepl("all", i, fixed = TRUE)){
    ss_all[[r]]=ss_all_old[[r]]
    }
  }
  if (!is.null(rsid) && rsid %in% colnames(ss_ov_old)){ss_ov$SNP=ss_ov_old[[rsid]]}
  if (!is.null(rsid) && rsid %in% colnames(ss_all_old)){ss_all$SNP=ss_all_old[[rsid]]}

  # Subset only intersection of chromosal position of each GWAS
  ss_ov$Allele1 = toupper(ss_ov$Allele1);ss_ov$Allele2 = toupper(ss_ov$Allele2);ss_all$Allele1 = toupper(ss_all$Allele1);ss_all$Allele2 = toupper(ss_all$Allele2)
  ss_tmp = inner_join(ss_ov,ss_all, by = c("CHR","POS","Allele1","Allele2"),suffix = c("_ov", "_all"))
  ss_tmp2 = inner_join(ss_ov,ss_all, by = c("CHR","POS","Allele1"="Allele2","Allele2"="Allele1"),suffix = c("_ov", "_all"))
  ss_tmp2 = ss_tmp2%>%mutate(BETA_all=(-1)*BETA_all,AF_Allele2_all=1-AF_Allele2_all)
  ss_tmp = rbind(ss_tmp,ss_tmp2)
  ss_tmp$P_all = as.numeric(ss_tmp$P_all) ; ss_tmp$P_all[ss_tmp$P_all==0] = 2.225074e-308
  ss_tmp$P_ov = as.numeric(ss_tmp$P_ov) ; ss_tmp$P_ov[ss_tmp$P_ov==0] = 2.225074e-308
  print("Number of variant of all, overlapped, subtracted samples")
  print(c(nrow(ss_all),nrow(ss_ov),nrow(ss_tmp)))

  # In case sample size is different
  for (i in c("n_all","n_ov","eff_all","eff_ov","n_case_all","n_ctrl_all","n_case_ov","n_ctrl_ov")){
    r=get(i)
    if (is.character(r)&grepl("ov", i, fixed = TRUE)){
    assign(i,ss_tmp[[r]])
    }
    if (is.character(r)&grepl("all", i, fixed = TRUE)){
    assign(i,ss_tmp[[r]])
    }
  }

  # Sample Size Calculation
  est_all = ss_tmp$BETA_all; se_all = ss_tmp$SE_all ;
  est_ov = ss_tmp$BETA_ov; se_ov =ss_tmp$SE_ov ;
  n_exc = n_all-n_ov

  # Get Minor Allele Frequency
  maf_all = ss_tmp$AF_Allele2_all
  maf_ov = ss_tmp$AF_Allele2_ov

  # Inverse Variance Weighting
  ss_tmp = ss_tmp%>%mutate(SE_IVW = sqrt(1/((1/SE_all^2)-(1/SE_ov^2))))
  ss_tmp = ss_tmp%>%mutate(BETA_IVW = (BETA_all/SE_all^2-BETA_ov/SE_ov^2)*SE_IVW^2)
  ss_tmp = ss_tmp%>%mutate(P_IVW = exp(pchisq((BETA_IVW^2)/SE_IVW^2,df=1,lower.tail=F,log.p=T)))

  # Reversed z-score
  if (!is.null(eff_all)&!is.null(eff_ov)) {
    w_all = sqrt(eff_all); w_ov = sqrt(eff_ov); w_RZ = sqrt(eff_all-eff_ov)
    z_all = PtoZ(ss_tmp$P_all,ss_tmp$BETA_all); z_ov = PtoZ(ss_tmp$P_ov,ss_tmp$BETA_ov)
    z_RZ = z_all*w_all - z_ov*w_ov ; z_RZ = z_RZ/w_RZ
    SE_RZ = 1/sqrt(2*maf_all*(1-maf_all)*((eff_all-eff_ov) + z_RZ^2))
    BETA_RZ =  z_RZ * SE_RZ
    sign_RZ = sign(z_RZ); P_RZ = exp(pnorm(abs(z_RZ),lower.tail=F,log.p=T))*2
  }
  # Partitioning GWAS
  if (phenotype=="binary"){
    xty_all = 0.5*n_all*maf_all-est_all*n_all*maf_all/4
    xty_ov = 0.5*n_ov*maf_ov-est_ov*n_ov*maf_ov/4
    xty_ex2  = xty_all - xty_ov
    est_exc2 = (-xty_ex2+0.5*(n_all*maf_all-n_ov*maf_ov))*4/(n_all*maf_all-n_ov*maf_ov)
    se_exc2 = sqrt((se_all^2*(n_all*maf_all)^2-se_ov^2*(n_ov*maf_ov)^2)/(n_exc*maf_all)^2)
  }
  if (phenotype=="continuous"){
    xty_all = est_all*n_all*maf_all
    xty_ov = est_ov*n_ov*maf_ov
    xty_ex2  = xty_all - xty_ov
    est_exc2 = xty_ex2/(n_all*maf_all-n_ov*maf_ov)
    se_exc2 = sqrt((se_all^2*(n_all*maf_all)^2-se_ov^2*(n_ov*maf_ov)^2)/(n_exc*maf_all)^2)
  }
  # Generate Output GWAS Dataframe
    new_ss = data.frame(
        CHR = ss_tmp$CHR, POS = ss_tmp$POS,
        Allele1 = ss_tmp$Allele1, Allele2 = ss_tmp$Allele2, AF_Allele = maf_all, N = n_exc, N_all = n_all, N_ov = n_ov,
        BETA_IVW = ss_tmp$BETA_IVW,SE_IVW = ss_tmp$SE_IVW,P_IVW = ss_tmp$P_IVW,
        BETA_Part = est_exc2,SE_Part = se_exc2,P_Part = exp(pchisq((est_exc2^2)/se_exc2^2,df=1,lower.tail=F,log.p=T)),
        BETA_all = ss_tmp$BETA_all,P_all = ss_tmp$P_all,SE_all = ss_tmp$SE_all,
        BETA_ov = ss_tmp$BETA_ov,P_ov = ss_tmp$P_ov, SE_ov = ss_tmp$SE_ov,
        xty_all = xty_all, xty_ov = xty_ov)
    if (!is.null(rsid)){new_ss$SNP = ss_tmp$SNP}
    if (!any(is.null(eff_all))&!any(is.null(eff_ov))){
        new_ss = new_ss %>% mutate(sign_RZ = sign_RZ, P_RZ = P_RZ, z_RZ = z_RZ, z_all = z_all, z_ov = z_ov, BETA_RZ = BETA_RZ, SE_RZ = SE_RZ, ESS_all = eff_all, ESS_ov = eff_ov, ESS_RZ = eff_all-eff_ov)

    }
  if (!dropna){
  new_ss$SE_IVW[is.na(new_ss$SE_IVW)] = max(ss_tmp$SE_all[is.na(new_ss$SE_IVW)],ss_tmp$SE_ov[is.na(new_ss$SE_IVW)])
  new_ss$SE_Part[is.na(new_ss$SE_Part)] = max(ss_tmp$SE_all[is.na(new_ss$SE_Part)],ss_tmp$SE_ov[is.na(new_ss$SE_Part)])
  }
  new_ss = new_ss %>% filter(!is.na(SE_IVW) & !is.na(SE_Part))
  # Check header and size of new overlap excluded GWAS
  print("New Summary")
  print(head(new_ss))
  print("Dimension of new summary")
  print(dim(new_ss))
  print(new_ss %>% filter(P_Part==0))
  # Write output
  fwrite(new_ss, output_file, col.names = T, row.names = F, quote = F,sep = " ",na = NA)
  #write.table(new_ss, output_file, col.names = T, row.names = F, quote = F)
  return(new_ss)
}

#' Overlap Adjustment
#'
#' Overlap Adjustment from a GWAS summary with summary of Overlapped samples
#' @param input_file Path of Summary statistics to be pruned
#' @param outdir Output directory of adjusted summary statistics
#' @param pt pvalue vector for pruning
#' @param pvals pvalue vector for Thresholding
#' @param hm_ref reference list
#' @return Output path of adjusted summary
#' @examples
#' Marker_select_ld(input_file='hm3_AS_0121.txt', outdir='./summary/', pt=c(1e-3,1e-4,1e-6,1e-8), pvals=c(0,0.01,0.05,0.1,0.2,0.5), hm_ref="eur")
#' @export
Marker_select_ld = function(input_file,outdir,pt=c(1e-3,1e-4,1e-6,1e-8),pvals=c(0,0.01,0.05,0.1,0.2,0.5),hm_ref){
if (grepl("eas",hm_ref,fixed = T)){
    ld=fread('/media/leelabsg-storage0/seokho/utils/Independent_LD/ldetect-data/ASN/fourier_ls-all.bed')
    ld = ld%>%mutate(CHR = as.numeric(substring(chr, 4)))
    ref_file = "/media/leelabsg-storage0/kisung/backup/kisung/GDA/reference/ldblk_1kg_eas/snpinfo_1kg_hm3"
}
if (grepl("eur",hm_ref,fixed = T)){
    ld=fread('/media/leelabsg-storage0/seokho/utils/Independent_LD/ldetect-data/EUR/fourier_ls-all.bed')
    ld = ld%>%mutate(CHR = as.numeric(substring(chr, 4)))
    ref_file = "/media/leelabsg-storage0/kisung/backup/kisung/GDA/reference/ldblk_1kg_eur/snpinfo_1kg_hm3"
}
if (grepl("eur",hm_ref,fixed = T)&grepl("hg38",hm_ref,fixed = T)){
    ld=fread('~/overlapGWAS/ADNI/EUR_LD_blocks.bed')
    ld = ld%>%mutate(CHR = as.numeric(substring(chr, 4)))
    ref_file = "/media/leelabsg-storage0/kisung/backup/kisung/GDA/reference/ldblk_1kg_eur/snpinfo_1kg_hm3"
}
if (grepl("eas",hm_ref,fixed = T)&grepl("hg38",hm_ref,fixed = T)){
    print("not yet built for eas-hg38")
    break
}

dir.create(outdir,showWarnings = F)
setwd(outdir)
ss_old = fread(input_file)
hm = fread(ref_file)
if(ishg38){
ss1 = inner_join(ss_old,hm,by = c("CHR","SNP","Allele1" = "A1","Allele2" = "A2"))
ss2 = inner_join(ss_old,hm,by = c("CHR","SNP","Allele1" = "A2","Allele2" = "A1"))
}else{
ss1 = inner_join(ss_old,hm,by = c("CHR","POS"="BP","Allele1" = "A1","Allele2" = "A2"))
ss2 = inner_join(ss_old,hm,by = c("CHR","POS"="BP","Allele1" = "A2","Allele2" = "A1"))
}
ss = rbind(ss1,ss2)
for (j in c("all","IVW","Part","RZ")){
    pname = paste0("P_",j); bname = paste0("BETA_",j)
    mkrs = (ss[ss[[pname]]<=pt,]) %>% select(CHR, POS)
    filt=NULL
    for (k in 1:nrow(mkrs)){
    filt = rbind(filt,ld %>% filter(CHR==as.numeric(mkrs[k,1]) & start<=as.numeric(mkrs[k,2]) & stop>=as.numeric(mkrs[k,2])))
    }
    filt = filt[!duplicated(filt),]
    sss = ss %>% select(CHR,POS) %>% mutate(CHRPOS = CHR*10^9+POS)
    filt = filt %>% mutate(chrstart = CHR*10^9+start,chrstop = CHR*10^9+stop) %>% select(chrstart,chrstop)
    sst = rep(NA,nrow(sss))
    for ( l in 1:nrow(sss)){
       sst[l] = !any(sss$CHRPOS[l] %[]% as.matrix(filt))
    }
    write.table(ss %>% select(SNP,A1=Allele1,A2=Allele2,BETA = !!bname,P=!!pname), paste0(basename(filename),'_',j,'_p9.txt'),col.names = T,row.names = F, quote = F)
    for ( i in pvals){
    tmp = ss[sst,]
    tmp = tmp[tmp[[pname]]>i & is.finite(tmp[[pname]]),]
    tmp = tmp %>% select(SNP,A1=Allele1,A2=Allele2,BETA = !!bname,P=!!pname)
    print(dim(tmp))
    write.table(tmp, paste0(basename(filename),'_',j,'_p',i,'.txt'),col.names = T,row.names = F, quote = F)
   }
}
if(!is.null(EraSOR)){
ss = ss %>% mutate(Chrpos=paste0(CHR,":",POS,"_",Allele1,"/",Allele2))
ss_er = fread(EraSOR)
ss = left_join(ss_er,ss, by = c("SNP"="Chrpos","CHR"))
mkrs = (ss[ss[["P"]]<=pt,]) %>% select(CHR, POS)
filt=NULL
for (i in 1:nrow(mkrs)){
filt = rbind(filt,ld %>% filter(CHR==as.numeric(mkrs[i,1]) & start<=as.numeric(mkrs[i,2]) & stop>=as.numeric(mkrs[i,2])))
}
filt = filt[!duplicated(filt),]
sss = ss %>% select(CHR,POS) %>% mutate(CHRPOS = CHR*10^9+POS)
filt = filt %>% mutate(chrstart = CHR*10^9+start,chrstop = CHR*10^9+stop) %>% select(chrstart,chrstop)
sst = rep(NA,nrow(sss))
for ( i in 1:nrow(sss)){
    sst[i] = !any(sss$CHRPOS[i] %[]% as.matrix(filt))
}
setwd(paste0("/home/seokhojeong/overlapGWAS/KoGES/1e",ldp,"/summary"))
write.table(ss %>% select(SNP=SNP.y,A1,A2,BETA=Z ,P),paste0(basename(filename),"_EraSOR_p9.txt"),col.names = T,row.names = F, quote = F)
for ( i in pvals){
ss_EraSOR = ss[sst,]
ss_EraSOR = ss_EraSOR %>% filter(P>i) %>% select(SNP=SNP.y,A1,A2,BETA=Z ,P)
print(dim(ss_EraSOR))
write.table(ss_EraSOR, paste0(basename(filename),'_EraSOR_p',i,'.txt'),col.names = T,row.names = F, quote = F)
}
}
}
# Sample usage
# source('/media/leelabsg-storage0/seokho/overlap/code/exclude_overlapped.r')

# When case/control sample size is not provided
# exclude_overlap("all.txt","ov.txt","out_n3.txt",n_all = 433540,n_ov = 72210,phenotype = "binary",dropna=F)
#############################################################################################################################################################
# Using case/control sample size and Effective sample size at the same time
# In this case, ESS of all samples are provided and case/control size of overlapped samples are provided
if(FALSE){
exclude_overlap("all_new.txt","/media/leelabsg-storage0/seokho/UKBB/overlapGWAS/Agen_T2D/KoGES_DM.tsv.gz","out_new_1004.txt",n_all = 433540,n_ov = 56918,
col_all = c("Beta","SE","EAF","Chr","Pos","Allele1","Allele2","P"),
eff_all =  211793, eff_ov = 10879.08,
col_ov = c("beta","sebeta","maf","chrom","pos","ref","alt","pval"),phenotype = "binary",dropna=T,
MAF = '/media/leelabsg-storage0/seokho/hapmap3/snpinfo_1kg_hm3')

exclude_overlap("T2D_ALL_Primary.txt","KBN_SAIGE_BINARY_DM.txt","hm3_0120.txt",n_all = 433540,n_ov = 56918,
col_all = c("Beta","SE","EAF","Chr","Pos","NEA","EA","P"),
col_ov = c("BETA","SE","AF_Allele2","chrom","POS","ref","alt","pval"),
eff_all =  211793, eff_ov = 10879.08,
phenotype = "binary",dropna=T,
MAF = NULL)

exclude_overlap("T2D_ALL_Primary.txt","/media/leelabsg-storage0/seokho/KoGES/DM/step2/DS_all_chr_locoT.txt","hm3_AS_0121.txt",n_all = 433540,n_ov = 5490,
col_all = c("Beta","SE","EAF","Chr","Pos","NEA","EA","P"),
eff_all =  211793, eff_ov =1279.278,
col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele2","Allele1","p.value"),phenotype = "binary",dropna=T,
MAF = NULL)

exclude_overlap("/media/leelabsg-storage0/seokho/UKBB/overlapGWAS/ADNI/IGAP_stage_1_mend.txt","/media/leelabsg-storage0/seokho/UKBB/overlapGWAS/ADNI/ADNI_data/summary_overlap_only_350.gz","subtracted.txt",n_all = 54162,n_ov = 327,
col_all = c("Beta","SE","NON","Chromosome","MarkerName","Non_Effect_allele","Effect_allele","pval"),
eff_all =  46668.53, eff_ov =279.2171,

col_ov = c("BETA","SE","AF_Allele2","CHR","MarkerID","Allele1","Allele2","p.value"),phenotype = "binary",dropna=T)

exclude_overlap("/home/seokhojeong/overlapGWAS/ADNI/Kunkle_etal_Stage1_results.txt","/home/seokhojeong/overlapGWAS/ADNI/summary_overlap_only_350.gz","/home/seokhojeong/overlapGWAS/ADNI/subtracted_new.txt",
n_all = 94437,n_ov = 327,
col_all = c("Beta","SE","NON","Chromosome","MarkerName","Non_Effect_allele","Effect_allele","Pvalue"),
eff_all = 88393.98, eff_ov =279.2171,
col_ov = c("BETA","SE","AF_Allele2","CHR","MarkerID","Allele1","Allele2","p.value"),phenotype = "binary",dropna=T,MAF = '/media/leelabsg-storage0/kisung/backup/kisung/GDA/reference/ldblk_1kg_eur/snpinfo_1kg_hm3')

exclude_overlap("all_new.txt","Korean_T2D_meta.txt","out_new_1220.txt",n_all = 433540,n_ov = 97676,
col_all = c("Beta","SE","EAF","Chr","Pos","Allele2","Allele1","P"),
eff_all =  211793, eff_ov = "n_samples",
col_ov = c("beta","se","eaf","chr","pos","reference_allele","other_allele","p-value"),phenotype = "binary",dropna=T,
MAF = '/media/leelabsg-storage0/seokho/hapmap3/snpinfo_1kg_hm3')

exclude_overlap("Kunkle_etal_Stage1_results_adni_common_grch38_hm3.txt", "summary_overlap_only_350.gz", "out_hm3.txt", col_all = c("Beta","SE","MAF","CHR","BP","Effect_allele","Non_Effect_allele","Pvalue"), col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele2","Allele1","p.value"), n_all = 63926, n_ov = 327, n_case_all = 21982, n_ctrl_all = 41944, n_case_ov = 226, n_ctrl_ov = 101, phenotype = "binary", dropna = T,rsid = "MarkerID")
}
# Korean_T2D_meta.txt.gz
#n_case_ov = 5083, n_ctrl_ov = 67127,
#exclude_overlap("ss_all.txt","ss_ov.txt","ss_out.txt",n_all = 326153,n_ov = 81538,
#n_case_all = 15156, n_ctrl_all =310996, n_case_ov = 3789, n_ctrl_ov = 77749,phenotype = "binary",dropna=F)
