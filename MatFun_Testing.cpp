#include <Rcpp.h>

using namespace Rcpp;
using namespace sugar;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericMatrix vecInsert(NumericVector v1, NumericVector v2, NumericVector v3) {
  NumericMatrix outmat;
  std::vector<double> tempVec;
  int len;

  if ( v1.length() != v2.length() || v1.length() != v3.length() ) {
    Rprintf("ERROR! Vector Lengths do not match!");
    R_FlushConsole();
    return outmat;
  }
  len = v1.length();
  for(int ii=0;ii<len;ii++){
    tempVec.push_back(rand());
  }
  outmat = NumericMatrix(4,v1.length());
  outmat(0,_)=Rcpp::as<NumericVector>(wrap(tempVec));

  return outmat;
}

// [[Rcpp::export]]
NumericVector upperVals(NumericVector v1) {
  int upper_id, v1_size;
  double upper;
  std::vector<double> outlist;
  
  v1_size = v1.size();
  upper_id = which_max(v1); //This is the index of the max element
  upper = floor(v1(upper_id));
  for (int ii=0;ii<v1_size;ii++) {
    if (v1(ii) >= upper) outlist.push_back(v1(ii));
  }
  return Rcpp::as<NumericVector>(wrap(outlist));
}

// [[Rcpp::export]]
NumericVector upper_ids(NumericVector v1) {
    int upper_id, v1_size;
  double upper;
  std::vector<int> outlist;
  
  v1_size = v1.size();
  upper_id = which_max(v1); //This is the index of the max element
  upper = floor(v1(upper_id));
  for (int ii=0;ii<v1_size;ii++) {
    if (v1(ii) >= upper) outlist.push_back(ii);
  }
  return Rcpp::as<NumericVector>(wrap(outlist));
}