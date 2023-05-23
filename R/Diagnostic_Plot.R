################################################################## Required Packages

##################################################################
#' Score Evaluation
#'
#' After applying polygenic risk score estimation tools, score evaluation for generating diagnostic plots
#' @importFrom data.table fread
#' @import ggplot2
#' @import dplyr
#' @param scores Evaluated personalized scores
#' @param title Plot title
#' @param Output_Plot_path Output Path for generated diagnostic plot
#' @param keep.ind Sample indices to be evaluated.
#' @param ref.ind Sample indices for Reference group
#' @param ref_name Reference group name
#' @param pheno_col Column name for phenotype (Optional)
#' @param covar_cols Specify PRS method platform (prscs, ldpred2)
#' @param method_names New labels for methods (Optional)
#' @return ggplot2 object of visual diagnostic
#' @examples
#' @export
diagnostic_plt <- function(scores,title,Output_Plot_path=NULL,keep.ind=NULL,ref.ind=NULL,ref_name = "",pheno_col=NULL,covar_cols=NULL,method_names=NULL){
  P<-auroc<-Method<-NULL
  if (is.character(scores)){
    ss_path = scores
    message("Loading input summary statistics")
    scores = fread(scores)
  }else if(!is.data.frame(scores)){
    stop("Invalid input format")
  }
  if(!is.null(keep.ind)){trgt_tmp=scores[keep.ind,]}else{trgt_tmp=scores}
  if(!is.null(ref.ind)){ref_tmp=scores[ref.ind,]
  ref = apply(ref_tmp %>% select(starts_with("prune")),2,function(x){as.numeric(auc(suppressMessages(roc(ref_tmp$pheno_label,x))))})
  }
  #if(!is.null(method_names)){method_names = paste0("Method ", nrow(nm))}
  trgt = apply(trgt_tmp %>% select(starts_with("prune")),2,function(x){as.numeric(auc(suppressMessages(roc(trgt_tmp$pheno_label,x))))})
  nm = t(as.data.frame(strsplit(names(trgt),split='_'),col.names = 1:length(trgt)))
  auc = data.frame(P=nm[,ncol(nm)],Method=rep(method_names,each=nrow(nm)/length(method_names)),auroc = trgt)
  auc$P[auc$P=='ref'] = -0.03
  auc$P = as.numeric(auc$P)
  auc$auroc = as.numeric(auc$auroc)
  auc_ref = auc
  if (!is.null(ref.ind)){
    auc_ref$auroc = as.numeric(ref)
    tmp = sort(unique(auc$P))
    auc$Method = factor(auc$Method,levels = method_names)
    p = ggplot(aes(x=P,y=auroc,col = Method),data = auc) + scale_color_brewer(palette = "Spectral") +
      geom_line(size = 1)+geom_point(size = 2.5)+labs(title =  title, x = "P-value", y = "AUROC")+ scale_x_continuous(breaks = tmp,labels = c("NP",as.character(tmp)[-1]),guide = guide_axis(check.overlap = TRUE))+ scale_y_continuous(breaks=c(0.45,0.5,0.6,0.7,0.8,0.9,1.0))+ theme_classic() +
      theme(plot.title = element_text(size = 35),legend.title = element_text(size=35),legend.text = element_text(size=30),axis.text = element_text(size = 20),axis.title = element_text(size = 30))
    p = p+ geom_line(aes(x=P,y=auroc,col = Method), linetype = 3,data = auc_ref,size = 0.75) + geom_text(aes(x = 0,y= 1,label=paste0("Pruned at P: ",nm[1,ncol(nm)-1])),col=1,check_overlap = T,size=6)
  }else{
    tmp = sort(unique(auc$P))
    auc$Method = factor(auc$Method,levels = method_names)
    p = ggplot(aes(x=P,y=auroc,col = Method),data = auc) + scale_color_brewer(palette = "Spectral") +
      geom_line(size = 1)+geom_point(size = 2.5)+labs(title = title, x = "P-value Threshold", y = "AUROC")+ scale_x_continuous(breaks = tmp,labels = c("NP",as.character(tmp)[-1]),guide = guide_axis(check.overlap = TRUE))+ scale_y_continuous(breaks=c(0.45,0.5,0.6,0.7,0.8,0.9,1.0))+ theme_classic() +
      theme(plot.title = element_text(size = 35),legend.title = element_text(size=35),legend.text = element_text(size=30),axis.text = element_text(size = 20),axis.title = element_text(size = 30)) + geom_text(aes(x = 0,y= 1,label=paste0("Pruned at P: ",nm[1,ncol(nm)-1])),col=1,check_overlap = T,size=6)
  }
  if(!is.null(Output_Plot_path)){
    ggsave(Output_Plot_path,p, width=18, height=4.5, units="in", dpi=450, pointsize=13)
  }
  print(p)
  return(p)
}
