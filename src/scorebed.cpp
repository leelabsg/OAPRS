#include <stdio.h>
#include <string>
#include <bitset>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Score Evaluation
//'
//' After applying polygenic risk score estimation tools, score evaluation for generating diagnostic plots
//'
//' @param fileName directory path of PRS result(e.g. prscs)
//' @param N Number of samples (fam)
//' @param P Number of variants (bim)
//' @param input Weight matrix for multiplication
//' @param extract logical vector indicating variants to be scored
//' @param keep logical vector indicating samples to be scored
//' @keywords internal
// [[Rcpp::export]]
arma::mat scorebed(const std::string fileName, int N, int P, const arma::mat input,
                   arma::Col<int> extract,  arma::Col<int> keep){
  std::ifstream bed;
  bed.open(fileName.c_str(), std::ios::binary);

  int m1; int m2; int m3;
  int colskip=0;
  int Nbytes=ceil(N/4);
  int n=(keep.n_elem>0) ? arma::sum(keep): N ;
  //int p=(extract.n_elem>0) ? arma::sum(extract): P ;
  char *header = new char[3];
  bed.read(header, 3);
  m1 = header[0] ; m2 = header[1]; m3 = header[2];
  try{
    if ((m1 != 108) | (m2 != 27)) throw 0;
  }catch(int exception)
  {
    bed.close();
    Rcpp::stop("Magic Number do not match! Not a plink format. \n");
  }
  printf("%d \n", m3);
  try{
    if (m3 != -1) throw 0;
  }catch(int exception)
  {
    bed.close();
    Rcpp::stop("Plink with individual-major mode is not supported \n");
  }
  std::bitset<8> bits;
  arma::mat scr = arma::mat(n, input.n_cols, arma::fill::zeros);
  arma::vec geno(n, arma::fill::zeros);
  char *char_blks = new char[Nbytes];
  int i=0;int j=0;int ii=0;int jj=0;
  float imp ; int nn; int gen=-1;
  while (i < P){
    if(extract[i]>0){
      bed.seekg(colskip,bed.cur);
      colskip=0;
      bed.read(char_blks, Nbytes);
      j=0;
      jj=0;
      imp=0; nn=0;
      while(j < Nbytes){
        bits = char_blks[j];
        //  std::string mystring =
        //    head.to_string<char,std::string::traits_type,std::string::allocator_type>();
        //  std::cout << "header "<< nn+1 << ": " << mystring << '\n';
        for (int l=0; l<4; l++){
          if (keep[j*4+l]>0){
            int g1 = bits[2*l];
            int g2 = bits[2*l+1];
            if ((g1==0) & (g2==0)){
              gen=2;imp+=gen;
            }
            if ((g1==1) & (g2==0)){
              gen=-1;geno[jj]=1;
              nn++;
            }
            if ((g1==1) & (g2==1)){
              gen=0;imp+=gen;
            }
            if ((g1==0) & (g2==1)){
              gen=1;imp+=gen;
            }
            for (unsigned int c=0; c<input.n_cols; c++){
              if(gen>=0){scr(jj,c) += input(ii,c)*gen;}
            }
            jj++;
          }
        }
        j++;
      }
     j=0;jj=0;imp = imp/(n-nn);
      while(j < Nbytes){
        for (int l=0; l<4; l++){
          if (keep[j*4+l]>0){
            if(geno[j*4+l]>0){
              for (unsigned int c=0; c<input.n_cols; c++){
              scr(jj,c) += input(ii,c)*imp;}
            }
            jj++;
          }
          }
        j++;
      }
      ii++;
    }
  else{
    colskip+=Nbytes;
  }
    i++;
  }
  bed.close();
  return scr;
}
