################################################################## Required Packages
globalVariables("..inds")
##################################################################

#' Score Evaluation
#'
#' After applying polygenic risk score estimation tools, score evaluation for generating diagnostic plots
#' @import data.table
#' @import dplyr
#' @param prs_res_paths directory path of PRS result(e.g. prscs)
#' @param lds Variant Filter results
#' @param target_path Plink binary file header of target/validation genotypes
#' @param output_score_path Path for storing output scores.(Optional)
#' @param platform Specify PRS method platform (prscs, ldpred2)
#' @param bim_path Optional target bim file path
#' @param fam_path Optional target fam file path
#' @param pheno_path Covariate file of target/validation samples for additional
#' @param pheno_col Column for phenotypes. If covariate file is not given. phenotype column from fam will be taken.
#' @param cl Cluster given from parallel
#' @return Formatted Summary Statistics
#' @examples
#' @export
score_eval <- function(prs_res_paths,lds,target_path,output_score_path=NULL,platform="prscs",bim_path=NULL,fam_path=NULL,pheno_path=NULL,pheno_col,cl=NULL){
  V6<-NULL
  methods=lds$methods
  pvals=lds$probs
  if(is.null(bim_path)){
    bim=fread(paste0(target_path,".bim"))
  }else{
    bim=fread(paste0(bim_path,".bim"))
  }
  if(is.null(fam_path)){
    fam=fread(paste0(target_path,".fam"))
  }else{
    fam=fread(paste0(fam_path,".fam"))
  }
  if(!is.null(pheno_path)){
    pheno = fread(pheno_path)
    fam = left_join(fam, pheno, by =c("V2"="DIST_ID")) %>% rename(pheno_label=!!pheno_col)
  }else{
    fam = fam %>% rename(pheno_label=V6)
  }
  scrs = data.frame(IID=fam$V2)
  for ( j in (1:length(methods))){
    method_name = methods[j]
    prs_res = fread(prs_res_paths[j])
    joint_ss = inner_join(prs_res,lds$selected_mkrs,by = "SNP")
    #
    prs_mkrs=inner_join(bim,joint_ss, by=c("V2"="SNP"),multiple="any")
    prs_mkrs$BETA_sign[(prs_mkrs$V5==prs_mkrs$A2 & prs_mkrs$V6==prs_mkrs$A1)] = -1
    prs_mkrs$BETA_sign[(prs_mkrs$V5==prs_mkrs$A1 & prs_mkrs$V6==prs_mkrs$A2)] = 1
    prs_mkrs$BETA_sign[!((prs_mkrs$V5==prs_mkrs$A1 & prs_mkrs$V6==prs_mkrs$A2)|(prs_mkrs$V5==prs_mkrs$A2 & prs_mkrs$V6==prs_mkrs$A1))] = 0

    prune_names = colnames(prs_mkrs)
    inds=grep(method_name,prune_names,fixed =T)
    weights=prs_mkrs$BETA_sign*prs_mkrs$BETA*prs_mkrs[,..inds]
    print(all(prs_mkrs$V2==bim$V2[bim$V2 %in% prs_mkrs$V2]))
    scr = apply(weights,2,pgs,bfile=target_path, extract=bim$V2 %in% prs_mkrs$V2,cluster=cl)
    scr = as.data.frame(scr) ; scr$FID=fam$V1; scr$IID=fam$V2
    scr = left_join(scr, fam, by = c("IID"="V2"));
    scrs = cbind(scrs,scr%>%select(starts_with("prune")))
  }
  scrs = left_join(scrs, fam, by = c("IID"="V2"));
  if(!is.null(output_score_path)){write.table(scrs,output_score_path,col.names = T, row.names = F, quote = F)}
  return(scrs)
}
