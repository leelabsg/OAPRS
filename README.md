# OAPRS

## ADNI adjusted IGAP summary 
#### [ADNI adjusted IGAP summary](https://storage.cloud.google.com/leelabsg/OAPRS/IGAP_summary_ADNI_overlap_adj.txt.gz?authuser=1)
You can download summary statistics generated using <code>OAPRS</code> on The International Genomics of Alzheimer's Project ([IGAP](https://www.niagads.org/datasets/ng00075)) \
with sample overlap adjustment of Alzheimer's Disease Neuroimaging Initiative ([ADNI](https://adni.loni.usc.edu/)) genotypes.

## Description
<code>OAPRS</code> is designed to help and guide adjusting sample overlap bias in building PRS without overfitting. \
<code>OAPRS</code> consists of four main steps: 1.summary information preparation, 2.sample overlap adjustment, 3.PRS construction, and 4.validation using visual diagnostics. \
Before the preparation step, GWAS summary statistics using only overlapped individual genotypes from the target data need to be generated using standard GWAS softwares. \
<code>OAPRS</code> assumes target data to be a PLINK binary file format.

## Detailed Manual
For list of functions and detailed options in <code>OAPRS</code>, Please refer to this manual : [OAPRS Manual](https://github.com/leelabsg/OAPRS/files/11571878/OAPRS.pdf)

## Dependencies
data.table, dplyr, ggplot2, RcppArmadillo, Rcpp (>= 1.0.9) 

R.utils are necessary for reading gzipped files <code> install.packages('R.utils') </code>

## Installation
You can install <code>OAPRS</code> Package from OAPRS github using [r-devtools](https://www.r-project.org/nosvn/pandoc/devtools.html)
```
devtools::install_github('leelabsg/OAPRS')
```

## Example Usage 
Example files are found in extdata of OAPRS package repository.

### 1. Format summary statistics 
Let's read a partial sample summary statistics on large scale genetic consortium with <code>Check_Sums</code>.
```
data_path = system.file("extdata/example",package = "OAPRS")
```
Create Column names first. By setting the cols variable, OAPRS reformat summary statistics for further use. 
```
cols = c(BETA="beta",Pval="pval",CHR="chrom",POS="pos",REF="ref",ALT="alt",SNP="rsids")
```
We can designate the summary statistics file path and genome_build, and population. If there is not specified column for sample size, Spcf_n, Spcf_n_case, Spcf_n_ctrl are needed.
```
cs = Check_Sums(paste0(data_path,'/consortium.ss'),
Genome_Build = "hg37", Pop = "eas",
cols=cols,
Spcf_n=249625,Spcf_n_case = 50466, Spcf_n_ctrl = 199159)
```
Similarly, we can format summary statistics with target summary.
```
cols = c(BETA="beta",Pval="pval",CHR="chrom",POS="pos",REF="ref",ALT="alt",SNP="rsids",SE="sebeta")
```
```
ts = Check_Sums(paste0(data_path,'/target.ss'),
Genome_Build = "hg37", Pop = "eas",
cols=cols,
Spcf_n=72210,Spcf_n_case = 5083, Spcf_n_ctrl = 62127)
```

### 2. Adjust consortium and target gwas 
With formatted summary statistics cs and ts, we can build adjusted summary statistics using <code>exclude_overlap</code>.

```
adj_ss = exclude_overlap(cs,ts,"adj.txt",phenotype="binary")
```
### 3. (External) Build PRS from prscs, lassosum, ldpred2, etc. 
In this example, <code>prscs</code> is applied In order to run prscs, we subset adjusted summary as input format.
```
library(dplyr)
write.table(adj_ss %>% select(SNP,A1=ALT,A2=REF,BETA=BETA_all,P=P_all),paste0(data_path,"/adj_ss_all.txt"),quote = F, col.names = T, row.names = F)
write.table(adj_ss %>% select(SNP,A1=ALT,A2=REF,BETA=BETA_IVW,P=P_IVW),paste0(data_path,"/adj_ss_IVW.txt"),quote = F, col.names = T, row.names = F)
write.table(adj_ss %>% select(SNP,A1=ALT,A2=REF,BETA=BETA_RZ,P=P_RZ),paste0(data_path,"/adj_ss_RZ.txt"),quote = F, col.names = T, row.names = F)
```
Here is some example usage for prscs for each adjustment methods.
```
for i in all IVW RZ
do
PRScs.py \
--ref_dir=snpinfo_1kg_hm3_eas.gz \
--bim_prefix=target \
--sst_file=adj_ss_${i}.txt \
--n_gwas=177415 \
--out_dir=${outdir}/${i}
done
```
### 4. Make filters for visual diagnostics 
Genome build and population information is essential in variant filtering for visual diagnostics.

```
lds = Marker_select_ld(adj_ss, Genome_Build="hg37", Pop = "eas")
```

### 5. Evaluate Scores on each thresholds
With result from 4 and prs weights on 3, scores will be evaluated by corresponding sample ID's.
You can specify platform for prs weights or specify column names 
```
prs_res_paths=c(paste0(data_path,'/prscs_all_chr_unadj.txt'),
  paste0(data_path,'/prscs_all_chr_IVW.txt'),paste0(data_path,'/prscs_all_chr_RZ.txt'))
scr = score_eval(prs_res_paths,lds, target_path = paste0(data_path,"/target"),platform="prscs",
  pheno_path = paste0(data_path,"/target.pheno"),ID_col="V2",pheno_col="V3")
```
For utilizing individual prs, columns with "ref" in the result can be used.
```
scr %>% select(IID, contains('ref'))
```
### 6. Generate Visual Diagnostics
Lets make visual diagnostics plot with generated scores. the legends are set using <code>method_names</code> option 
```
diagnostic_plt(score=scr,title="test",
  Output_Plot_path="~/plot.png",method_names = c("all","IVW","RZ"))
```

