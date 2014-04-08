#include <Rcpp.h>
#include <queue>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector priorityQueueSort(NumericVector invals) {
   NumericVector outvals(invals.size());
   std::priority_queue<double> sortQueue;
   int i;
   i=0;
   for(NumericVector::iterator it=invals.begin(); it != invals.end(); ++it) {
     sortQueue.push(*it);
   }
   while(!sortQueue.empty()) {
     outvals[i]=sortQueue.top();
     sortQueue.pop();
     i++;
   }
   return outvals;
}
