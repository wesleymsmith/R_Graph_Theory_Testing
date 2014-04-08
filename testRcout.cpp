#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
void testRcout(NumericVector vl) {
   double outnum;
   for (int ii=0; ii<vl.length();ii++) {
     outnum=vl[ii];
     Rprintf("entery at index %d is %2.1f \n ",ii,outnum);
     R_FlushConsole();
   }
   return;
}
