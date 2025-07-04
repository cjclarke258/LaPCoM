#include <Rcpp.h>
using namespace Rcpp;

// Function to compute factorial
double factorial(int n) {
  if (n == 0)
    return 1;
  else
    return n * factorial(n - 1);
} // end factorial function

// Main function
// [[Rcpp::export]]
double log_LPM_fast(NumericMatrix Y_m, double alpha, NumericMatrix D, std::string net_type, std::string net_mode) {
  
  // needed no matter the value of net_type
  int N = Y_m.nrow();
  
  // initialise the log-likelihood return value
  double ll = 0;

  // net_mode options
  if (net_mode == "undirected") {
    // net_type options
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
          if (j != i) {
            if (net_type == "binary") {
              ll += ( (alpha - D(i, j)) * (Y_m(i, j)) ) - ( (log(1 + exp(alpha - D(i, j)))) ); // updated
            } else if (net_type ==  "count") {
              // ll += (( Y_m(i, j) * (alpha - D(i, j)) ) - ( exp(alpha - D(i, j)) ) - ( log( factorial(Y_m(i, j)) ) ));
              ll += (( Y_m(i, j) * (alpha - D(i, j)) ) - ( exp(alpha - D(i, j)) ) );
            } else if (net_type == "pos_real") {
              ll += ( ( log(- (1 / (alpha - D(i, j)))) ) + (Y_m(i, j) / (alpha - D(i, j))) );
            } else {
              Rcpp::Rcout << "Uh oh. We only support binary, count or pos_real networks right now." << std::endl;
            } // end if else statement
          } // end if statement
        } // end j for loop
      } // end i for loop
  } else if (net_mode == "directed") {
      // net_type options
      for (int i = 0; i < N; i++) {
          for (int j = 0; j < N; j++) {
            if (j != i) {
              if (net_type == "binary") {
                ll += ( (alpha - D(i, j)) * (Y_m(i, j)) ) - ( (log(1 + exp(alpha - D(i, j)))) ); // updated
              } else if (net_type ==  "count") {
                // ll += (( Y_m(i, j) * (alpha - D(i, j)) ) - ( exp(alpha - D(i, j)) ) - ( log( factorial(Y_m(i, j)) ) ));
                ll += (( Y_m(i, j) * (alpha - D(i, j)) ) - ( exp(alpha - D(i, j)) ) );
              } else if (net_type == "pos_real") {
                ll += ( ( log(- (1 / (alpha - D(i, j)))) ) + (Y_m(i, j) / (alpha - D(i, j))) );
              } else {
                Rcpp::Rcout << "Uh oh. We only support binary, count or pos_real networks right now." << std::endl;
              } // end if else statement
            } // end if statement
          } // end j for loop
        } // end i for loop
  } else {
      Rcpp::Rcout << "Uh oh. You must indicate whether your multiplex contains 'undirected' or 'directed' networks." << std::endl;
  } // end if else statement
  
  return ll;
} // end log_LPM_cpp function