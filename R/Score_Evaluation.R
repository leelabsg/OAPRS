################################################################## Required Packages
globalVariables("..inds")
##################################################################

#' Score Evaluation
#'
#' After applying polygenic risk score estimation tools, score evaluation for generating diagnostic plots
#' @import data.table
#' @import dplyr
#' @importFrom lassosum pgs
#' @param prs_res_paths directory path of PRS result(e.g. prscs)
#' @param lds Variant Filter results
#' @param target_path Plink binary file header of target/validation genotypes
#' @param output_score_path Path for storing output scores.(Optional)
#' @param platform Specify PRS method platform (prscs, ldpred2, PRsice2)
#' @param prs_cols For other platforms, you can provide colnames with "SNP","POS","A1","A2","BETA".
#' @param bim_path Optional target bim file path
#' @param fam_path Optional target fam file path
#' @param pheno_path Covariate file of target/validation samples for additional
#' @param ID_col Column name for Sample ID in pheno_file
#' @param pheno_col Column name for phenotypes. If covariate file is not given. phenotype column from fam will be taken.
#' @return Formatted Summary Statistics
#' @examples
#' @export
score_eval <- function(prs_res_paths,lds,target_path,output_score_path=NULL,platform="prscs",prs_cols=NULL,bim_path=NULL,fam_path=NULL,pheno_path=NULL,ID_col,pheno_col){
  V6<-NULL
  methods=lds$methods
  pvals=lds$probs
  if(is.null(bim_path)){
    bim=fread(paste0(target_path,".bim"))
  }else{
    bim=fread(paste0(bim_path,".bim"))
  }
  P=nrow(bim)
  if(is.null(fam_path)){
    fam=fread(paste0(target_path,".fam"))
  }else{
    fam=fread(paste0(fam_path,".fam"))
  }
  N=nrow(fam)
  if(!is.null(pheno_path)){
    pheno = fread(pheno_path)
    fam = left_join(fam, pheno, by =c("V2"=ID_col),suffix = c(".fam", "")) %>% rename(pheno_label=!!pheno_col)
  }else{
    fam = fam %>% rename(pheno_label=V6)
  }
  weights = NULL
  exc_mkrs = c()
  for ( j in (1:length(methods))){
    method_name = methods[j]
    prs_res = fread(prs_res_paths[j])
    if (platform%in%c("prscs","ldpred2")){colnames(prs_res) = c("CHR","SNP","POS","A1","A2","BETA")
    }else if (platform=="PRsice2 "){
      colnames(prs_res) = c("CHR","SNP","POS","A1","A2","BETA","P")
    }else if (is.null(prs_cols)|length(prs_cols)!=ncol(prs_res)){
      stop("Please specify columns for prs file")
    }else{colnames(prs_res) = prs_cols}
    joint_ss = inner_join(prs_res,lds$selected_mkrs,by = c("SNP"))
    prs_mkrs=inner_join(bim,joint_ss, by=c("V2"="SNP"),multiple="any")
    prs_mkrs$BETA_sign[(prs_mkrs$V5==prs_mkrs$A2 & prs_mkrs$V6==prs_mkrs$A1)] = -1
    prs_mkrs$BETA_sign[(prs_mkrs$V5==prs_mkrs$A1 & prs_mkrs$V6==prs_mkrs$A2)] = 1
    prs_mkrs$BETA_sign[!((prs_mkrs$V5==prs_mkrs$A1 & prs_mkrs$V6==prs_mkrs$A2)|(prs_mkrs$V5==prs_mkrs$A2 & prs_mkrs$V6==prs_mkrs$A1))] = 0
    prune_names = colnames(prs_mkrs)
    inds=grep(method_name,prune_names,fixed =T)
    weights=cbind(weights,prs_mkrs$BETA_sign*prs_mkrs$BETA*prs_mkrs[,..inds])
    exc_mkrs = c(exc_mkrs,prs_mkrs$V2)
  }
  exc_mkrs = unique(exc_mkrs)
    scr = scorebed(fileName=paste0(target_path,".bed"), N=N,P=P,extract=as.numeric(bim$V2 %in% exc_mkrs),input = as.matrix(weights),keep=rep(1,N))
    colnames(scr) = colnames(weights)
    scr = as.data.frame(scr) ; scr$FID=fam$V1; scr$IID=fam$V2
    scr = left_join(scr, fam, by = c("IID"="V2"));
  if(!is.null(output_score_path)){write.table(scr,output_score_path,col.names = T, row.names = F, quote = F)}
  return(scr)
}
