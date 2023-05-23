# OAPRS
Package for sample Overlap adjustment in PRS
## ADNI adjusted IGAP summary 
[download](https://storage.cloud.google.com/leelabsg/OAPRS/IGAP_summary_ADNI_overlap_adj.txt.gz?authuser=1)

## Dependencies
Lassosum \
R.utils \
data.table \
dplyr \
Desctools 

## Installation
To install lassosum,
```
devtools::install_github('tshmak/lassosum')
```
Install OAPRS Package 
```
devtools::install_github('leelabsg/OAPRS')
```

## Example Usage 

### 1. Format summary statistics 
Let's read a sample summary statistics on large scale genetic consortium.
```
ss_all=system.file('extdata/example/consortium.ss',package = "OAPRS")

cs = Check_Sums(system.file('extdata/example/consortium.ss',package = 'OAPRS'),
Genome_Build = "hg37", Pop = "eas",
cols = c(BETA="beta",Pval="pval",CHR="chrom",POS="pos",REF="ref",ALT="alt",SNP="rsids"),
Spcf_n=249625,Spcf_n_case = 50466, Spcf_n_ctrl = 199159)
```
Similarly, we can format summary statistics with target summary.
```
ss_trgt=system.file('extdata/example/target.ss',package = "OAPRS")
ts = Check_Sums(system.file('extdata/example/target.ss',package = 'OAPRS'),
Genome_Build = "hg37", Pop = "eas",
cols = c(BETA="beta",Pval="pval",CHR="chrom",POS="pos",REF="ref",ALT="alt",SNP="rsids",SE="sebeta"),
Spcf_n=72210,Spcf_n_case = 5083, Spcf_n_ctrl = 62127)
```

### 2. Adjust consortium and target gwas 
```
adj_ss = exclude_overlap(cs,ts,"adj.txt",phenotype="binary")
```
### 3. (External) Build PRS from prscs, lassosum, ldpred2, etc. 
Format for prscs
```
write.table(adj_ss %>% select(SNP,A1=ALT,A2=REF,BETA=BETA_all,P=P_all),"adj_ss_all.txt",quote = F, col.names = T, row.names = F)
write.table(adj_ss %>% select(SNP,A1=ALT,A2=REF,BETA=BETA_IVW,P=P_IVW),"adj_ss_IVW.txt",quote = F, col.names = T, row.names = F)
write.table(adj_ss %>% select(SNP,A1=ALT,A2=REF,BETA=BETA_RZ,P=P_RZ),"adj_ss_RZ.txt",quote = F, col.names = T, row.names = F)
```

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
```
lds = Marker_select_ld(adj_ss, Genome_Build="hg37", Pop = "eas")
```

### 5. Evaluate Scores on each thresholds
```

```

### 6. Generate Visual Diagnostics
