#include <Rcpp.h>

using namespace Rcpp;

//Function unique = Environment::base_env["unique"];

double biggo = std::numeric_limits<double>::max();

double max(double x1, double x2) {
  if (x1 > x2) return(x1);
  return(x2);
}


// [[Rcpp::export(.dmat)]]

NumericMatrix dmat(NumericVector x1, NumericVector x2, NumericVector y1, NumericVector y2) {
  int n1 = x1.size();
  int n2 = x2.size();
  if (y1.size() != n1) Rcpp::stop("X and Y lengths differ for row points.");
  if (y2.size() != n2) Rcpp::stop("X and Y lengths differ for column points.");
  NumericMatrix dists(n1,n2); 
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      dists(i1,i2) = sqrt(pow(x1[i1] - x2[i2],2) + pow(y1[i1] - y2[i2],2));
    }
  }
  return(dists);
}            

// [[Rcpp::export(.dmatex)]]

NumericMatrix dmatex(NumericVector x1, NumericVector x2, NumericVector y1, NumericVector y2, double pwr) {
  int n1 = x1.size();
  int n2 = x2.size();
  if (y1.size() != n1) Rcpp::stop("X and Y lengths differ for row points.");
  if (y2.size() != n2) Rcpp::stop("X and Y lengths differ for column points.");
  NumericMatrix dists(n1,n2); 
  if (pwr == R_PosInf) {
    for (int i1 = 0; i1 < n1; ++i1) {
      for (int i2 = 0; i2 < n2; ++i2) {
        dists(i1,i2) = max(std::abs(x1[i1] - x2[i2]),std::abs(y1[i1] - y2[i2]));
      }
    }
    return(dists);
  }
  for (int i1 = 0; i1 < n1; ++i1) {
    for (int i2 = 0; i2 < n2; ++i2) {
      dists(i1,i2) = pow(pow(std::abs(x1[i1] - x2[i2]),pwr) + pow(std::abs(y1[i1] - y2[i2]),pwr),1.0/pwr);
    }
  }
  return(dists);
} 


// [[Rcpp::export(.rviss)]]

IntegerVector rviss(NumericMatrix dm, IntegerVector ss) {
  int n1 = dm.nrow(); // Number of rows
  int n3 = dm.ncol();
  int n2 = ss.size(); // Number of cols in subset
  IntegerVector result(n1); // Loop through matrix,  finding smallest
  for (int i = 0; i < n1; ++i) {
    double dmin = biggo;
    int imin = 0;
    for (int j = 0; j < n2; ++j) {
      if (ss[j] > n3) Rcpp::stop("Index of location out of range.");
      if (dm(i,ss[j] - 1) < dmin) {
        dmin = dm(i,ss[j] - 1);
        imin = clone(ss)[j];
      }
    }
    result[i] = imin;  // Correct for R index starting at 1, but cpp at 0
  }
  return result;
}


// [[Rcpp::export(.dtotal)]]

double dtotal(NumericMatrix dm, IntegerVector ss) {
  int n1 = dm.nrow(); // Number of rows
  IntegerVector nni(n1);
  nni = rviss(dm, ss);
  double sum = 0.0;
  for (int i = 0; i < n1; ++i) {
    sum += dm(i,nni[i] - 1);
  }
  return sum;
}


// [[Rcpp::export(.bestswap)]]

IntegerVector bestswap(NumericMatrix dm, IntegerVector ins, IntegerVector outs) {
  int n1 = ins.size();
  int n2 = outs.size();
  IntegerVector testins = clone(ins);
  double dbest = dtotal(dm,testins);
  IntegerVector result = clone(testins);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
       testins = clone(ins);
       testins[i] = outs[j];
       double dtest = dtotal(dm,testins);
       if (dtest < dbest) {
         dbest = dtest;
         result = testins;
       }
    }
  }
  return result;
}

// [[Rcpp::export(.bestswap2)]]

IntegerVector bestswap2(NumericMatrix dm, IntegerVector ins, IntegerVector outs, int n_force) {
  int n1 = ins.size();
  int n2 = outs.size();
  IntegerVector testins = clone(ins);
  double dbest = dtotal(dm,testins);
  IntegerVector result = clone(testins);
  for (int i = n_force; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
       testins = clone(ins);
       testins[i] = outs[j];
       double dtest = dtotal(dm,testins);
       if (dtest < dbest) {
         dbest = dtest;
         result = testins;
       }
      }
    }
  return result;
}

// Might get used later

int p1med(NumericMatrix dm) {
  int n1 = dm.nrow();
  int n2 = dm.ncol();
  double best = biggo;
  int bestIndex = 0;
  for (int i = 0; i < n2; ++i) {
    double rowsum = 0.0;
    for (int j = 0; j < n1; ++j) rowsum += dm(j,i);
    if (rowsum < best) {
      best = rowsum;
      bestIndex = i + 1;
    }
  }
  return bestIndex;
} 




