################################################################## Required Packages

##################################################################
#' Variant filtering for generating diagnostic plots
#'
#' Variant filtering by Independent LD blocks
#' @import data.table
#' @import dplyr
#' @param ss Summary statistic object or path for formatted summary statistics
#' @param pt Threshold for pruning
#' @param probs Multiple P-value cutoffs for diagnostic plot
#' @param Genome_Build Genome build. hg37 or hg38
#' @param Pop Population group (EUR or EAS)
#' @param methods Specify sample size if not stated in column
#' @return Dataframe of logical filter of markers, threshold, methods
#' @examples
#' @export
Marker_select_ld = function(ss,pt=1e-4,probs=c(0,0.01,0.05,0.1,0.2,0.5),Genome_Build,Pop,methods= c("all","IVW","RZ")){
  MAF<-CHR<-POS<-chr<-SNP<-NULL
  if (is.character(ss)){
    ss_path = ss
    message("Loading input summary statistics")
    ss = fread(ss)
  }else if(!is.data.frame(ss)){
    stop("Invalid input format")
  }

  if (grepl("eas",Pop,fixed = T)&grepl("hg37",Genome_Build,fixed = T)){
    ld=fread(system.file('extdata/GRCh37_blks/fourier_ls_ASN.bed',package="OAPRS"))
    ld = ld%>%mutate(CHR = as.numeric(substring(chr, 4)))
    ref_file = system.file("extdata/SNPInfo/snpinfo_1kg_hm3_eas.gz",package="OAPRS")
  }
  if (grepl("eur",Pop,fixed = T)&grepl("hg37",Genome_Build,fixed = T)){
    ld=fread(system.file('extdata/GRCh37_blks/fourier_ls_EUR.bed',package="OAPRS"))
    ld = ld%>%mutate(CHR = as.numeric(substring(chr, 4)))
    ref_file = system.file("extdata/SNPInfo/snpinfo_1kg_hm3_eur.gz",package="OAPRS")
  }
  if (grepl("eur",Pop,fixed = T)&grepl("hg38",Genome_Build,fixed = T)){
    ld=fread(system.file('extdata/GRCh38_blks/deCODE_EUR_LD_blocks.bed',package="OAPRS"))
    ld = ld%>%mutate(CHR = as.numeric(substring(chr, 4)))
    ref_file = system.file("extdata/SNPInfo/snpinfo_1kg_hm3_eur.gz",package="OAPRS")
  }
  if (grepl("eas",Pop,fixed = T)&grepl("hg38",Genome_Build,fixed = T)){
    ld=fread(system.file('extdata/GRCh38_blks/pyrho_EAS_LD_blocks.bed',package="OAPRS"))
    ld = ld%>%mutate(CHR = as.numeric(substring(chr, 4)))
    ref_file = system.file("extdata/SNPInfo/snpinfo_1kg_hm3_eas.gz",package="OAPRS")
  }

  mkr_select_int = function(ss,pname, pt,ld){
    by_chr = function(chrom,ss,ld){
      ld_tmp=ld %>% filter(CHR==chrom)
      cut_ls=c(ld_tmp$start,ld_tmp$stop[nrow(ld_tmp)])
      ss_tmp=ss %>% filter(CHR==chrom) %>% mutate(blk_lab=cut(POS,breaks=cut_ls, label=F))
      mkrs = unique(ss_tmp$blk_lab[ss_tmp[[pname]]<=pt])
      return(ss_tmp$CHR_POS[!ss_tmp$blk_lab%in%mkrs])
    }
    filtered_list=unlist(sapply(1:22,by_chr,ld=ld,ss=ss))
    sst = ss$CHR_POS %in% filtered_list
    return(sst)
  }
  hm = fread(ref_file)

  ss1 = inner_join(ss,hm%>%select(-MAF),by = c("CHR","POS"="BP","REF" = "A1","ALT" = "A2","SNP"))
  ss2 = inner_join(ss,hm%>%select(-MAF),by = c("CHR","POS"="BP","REF" = "A2","ALT" = "A1","SNP"))
  ss = rbind(ss1,ss2)

  selected_mkrs = ss %>% select(CHR,POS,SNP)
  sst = data.frame()
  for (j in methods){
    pname = paste0("P_",j)
    filt=mkr_select_int(ss,pname, pt, ld)
    sst=setDT(sst)[, (paste0("prune_",j,"_",pt,"_ref")) := ss[[pname]]>0]
    for (prob in probs){
      sst=setDT(sst)[, (paste0("prune_",j,"_",pt,"_",prob)) := (filt & ss[[pname]]>prob)]
    }
  }
  lds = list(selected_mkrs=cbind(selected_mkrs,sst),prune = pt,probs = probs,methods=methods)
  print(colSums(lds$selected_mkrs[,4:ncol(lds$selected_mkrs)]))
  return(lds)
}

