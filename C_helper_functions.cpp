#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List pc_sim(NumericVector pars, int its) {
  
  
  NumericVector range = {0};
  NumericVector qrange = Rcpp::pnorm(range, pars[0], pars[1]);
  NumericVector rrunif = Rcpp::runif(its, qrange[0], 1);
  NumericVector rspeed = Rcpp::qnorm(rrunif, pars[0], pars[1]);
  
  NumericVector nreads = Rcpp::rpois(its, 1000);
  //NumericVector rspeed = Rcpp::rnorm(its, pars[0], pars[1]);
  NumericVector superdist (sum((nreads * nreads - nreads) / 2) / 2);
  int K = 0;
  
  for(int it = 0; it < its; ++it) {
    
    NumericVector reads  = Rcpp::runif(nreads[it], -200, 200);
    LogicalVector reads2 = ((reads >= 30*rspeed[it]) & (reads <= 45*rspeed[it])) | ((reads <= -30*rspeed[it]) & (reads >= -45*rspeed[it])) ;
    if(is_true(any(is_na(reads2)))) break;
    NumericVector x = reads[reads2];
    
    int n = x.size();
    
    if(n < 2) break;
    
    NumericVector from ((n * n - n) / 2);
    NumericVector to   ((n * n - n) / 2);
    int k = 0;
    
    for(int i = 0; i < n; ++i) {
      for(int j = i + 1; j < n; ++j){
        
        from[k] = x[i];
        to[k]   = x[j];
        
        ++k;
      }
    }
    
    NumericVector dist   = abs(to - from);
    LogicalVector filter = dist > (15 * rspeed[it]);
    if(is_true(any(is_na(filter)))) break;
    NumericVector dist2  = dist[filter];
    
    K = K + dist2.size();
    
    if(dist2.size() == 0) {
      break;
    } else if(dist2.size() == 1) {
      superdist[K - 1] = dist2[0];
    } else {
      superdist[Rcpp::Range(K - dist2.size(), K - 1)] = dist2;
    }
  }
  
  if(K == 0) {
    List out = List::create(Named("u") = 0, Named("s") = 0);
    return out;
  }
  
  double uest = Rcpp::mean(superdist[Rcpp::Range(0, K-1)]) / 75;
  double sest = Rcpp::sd(superdist[Rcpp::Range(0, K-1)]) / 75;
  
  List out = List::create(Named("u") = uest, Named("s") = sest);
  
  return out;
  
}

// [[Rcpp::export]]
List dist_1d(IntegerVector x, int MAX) {
  
  int n = x.size();
  IntegerVector from ((n * n - n) / 2);
  IntegerVector to   ((n * n - n) / 2);
  int k = 0;
  
  for(int i = 0; i < n; ++i) {
    for(int j = i + 1; j < n; ++j){
      
      from[k] = x[i];
      to[k]   = x[j];
      
      ++k;
    }
  }
  
  IntegerVector dist = to - from; 
  LogicalVector filter = dist < MAX;
  IntegerVector from_f = from[filter];
  IntegerVector dist_f = dist[filter];
  
  List final =  List::create(Named("from") = from_f, Named("dist") = dist_f);
  
  return final;
}

// [[Rcpp::export]]
List EM_doublemix(NumericVector x, NumericVector parprior, double tol = 1e-6, double minexppar = 5000) {
  
  int n = x.size();
  NumericVector EM_exp (n);
  NumericVector EM_unif_v (n);
  NumericVector EM_normalise (n);
  NumericVector EM_parprior = parprior;
  NumericVector all_loglik (100);
  double EM_unif_d   = 0;
  double EM_exp_sum  = 0;
  double EM_unif_sum = 0;
  
  // first E step
  
  EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[2];
  EM_unif_d     =  1 / (parprior[1]) * parprior[3];
  
  EM_normalise  = EM_exp + EM_unif_d;
  
  EM_exp        = EM_exp  / EM_normalise;
  EM_unif_v     = EM_unif_d / EM_normalise;
  
  all_loglik[0] = sum(log(EM_normalise));
  
  for(int i = 1; i < 100; ++i) {
    
    // M step, update parameters
    EM_exp_sum  = sum(EM_exp);
    EM_unif_sum = sum(EM_unif_v);
    
    EM_parprior[0] = sum( x * EM_exp) / EM_exp_sum;
    
    // M step, update priors
    EM_parprior[2] = EM_exp_sum  / n;
    EM_parprior[3] = EM_unif_sum / n;
    
    // M step, restrict parameters
    if (EM_parprior[0] < minexppar) EM_parprior[0] = minexppar;
    if (EM_parprior[2] < 0.03) EM_parprior[2] = 0.03;
    if (EM_parprior[3] > 0.97) EM_parprior[3] = 0.97;
    
    // E step, update posterior
    EM_exp    = (1/EM_parprior[0] * exp(-1/EM_parprior[0] * x)) * EM_parprior[2];
    EM_unif_d =  1 / (EM_parprior[1]) * EM_parprior[3];
    
    EM_normalise = EM_exp + EM_unif_d;
    
    EM_exp    = EM_exp  / EM_normalise;
    EM_unif_v = EM_unif_d / EM_normalise;
    
    all_loglik[i] = sum(log(EM_normalise));
    
    if ((all_loglik[i] - all_loglik[i - 1]) / -all_loglik[i - 1] < tol) break;
    
  }
  
  List final = List::create(Named("pars_priors") = EM_parprior, Named("loglike") = all_loglik);
  
  return  final;    
}   

// [[Rcpp::export]]
List EM_triplemix(NumericVector x, NumericVector parprior, double tol = 1e-6, double minexppar = 5000) {
  
  int n = x.size();
  NumericVector EM_exp (n);
  NumericVector EM_norm (n);
  NumericVector EM_unif_v (n);
  NumericVector EM_normalise (n);
  NumericVector EM_parprior = parprior;
  NumericVector all_loglik (100);
  double EM_unif_d   = 0;
  double EM_exp_sum  = 0;
  double EM_norm_sum = 0;
  double EM_unif_sum = 0;
  
  // first E step
  
  EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[4];
  EM_norm   = exp(-pow(x - parprior[1], 2) / 2 / parprior[2]) / (2.5066282746 * sqrt(parprior[2])) * parprior[5];
  EM_unif_d     =  1 / (parprior[3]) * parprior[6];
  
  EM_normalise  = EM_exp + EM_norm + EM_unif_d;
  
  EM_exp        = EM_exp  / EM_normalise;
  EM_norm       = EM_norm / EM_normalise;
  EM_unif_v     = EM_unif_d / EM_normalise;
  
  all_loglik[0] = sum(log(EM_normalise));
  
  for(int i = 1; i < 100; ++i) {
    
    // M step, update parameters
    EM_exp_sum  = sum(EM_exp);
    EM_norm_sum = sum(EM_norm);
    EM_unif_sum = sum(EM_unif_v);
    
    EM_parprior[0] = sum( x * EM_exp) / EM_exp_sum;
    EM_parprior[1] = sum( x * EM_norm) / EM_norm_sum;
    EM_parprior[2] = sum(pow(x - EM_parprior[1], 2) * EM_norm) / EM_norm_sum;
    
    // M step, update priors
    EM_parprior[4] = EM_exp_sum  / n;
    EM_parprior[5] = EM_norm_sum / n;
    EM_parprior[6] = EM_unif_sum / n;
    
    // M step, restrict parameters
    if (EM_parprior[0] < minexppar) EM_parprior[0] = minexppar;
    if (EM_parprior[4] < 0.05) EM_parprior[4] = 0.05;
    if (EM_parprior[6] > 0.95) EM_parprior[6] = 0.95;
    
    // E step, update posterior
    EM_exp    = (1/EM_parprior[0] * exp(-1/EM_parprior[0] * x)) * EM_parprior[4];
    EM_norm   = exp(-pow(x - EM_parprior[1], 2) / 2 / EM_parprior[2]) / (2.5066282746 * sqrt(EM_parprior[2])) * EM_parprior[5];
    EM_unif_d =  1 / (EM_parprior[3]) * EM_parprior[6];
    
    EM_normalise = EM_exp + EM_norm + EM_unif_d;
    
    EM_exp    = EM_exp  / EM_normalise;
    EM_norm   = EM_norm / EM_normalise;
    EM_unif_v = EM_unif_d / EM_normalise;
    
    all_loglik[i] = sum(log(EM_normalise));
    
    if ((all_loglik[i] - all_loglik[i - 1]) / -all_loglik[i - 1] < tol) break;
    
  }
  
  List final = List::create(Named("pars_priors") = EM_parprior, Named("loglike") = all_loglik);
  
  return  final;    
}   

// [[Rcpp::export]]
List EM_quadruplemix(NumericVector x, NumericVector parprior, double tol = 1e-6, double minexppar = 5000) {
  
  int n = x.size();
  NumericVector EM_exp (n);
  NumericVector EM_norm1 (n);
  NumericVector EM_norm2 (n);
  NumericVector EM_unif_v (n);
  NumericVector EM_normalise (n);
  NumericVector EM_parprior = parprior;
  NumericVector all_loglik (100);
  double EM_unif_d   = 0;
  double EM_exp_sum  = 0;
  double EM_norm1_sum = 0;
  double EM_norm2_sum = 0;
  double EM_unif_sum = 0;
  
  // first E step
  
  EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[6];
  EM_norm1      = exp(-pow(x - parprior[1], 2) / 2 / parprior[2]) / (2.5066282746 * sqrt(parprior[2])) * parprior[7];
  EM_norm2      = exp(-pow(x - parprior[3], 2) / 2 / parprior[4]) / (2.5066282746 * sqrt(parprior[4])) * parprior[8];
  EM_unif_d     =  1 / (parprior[5]) * parprior[9];
  
  EM_normalise  = EM_exp + EM_norm1 + EM_norm2 + EM_unif_d;
  
  EM_exp        = EM_exp  / EM_normalise;
  EM_norm1      = EM_norm1 / EM_normalise;
  EM_norm2      = EM_norm2 / EM_normalise;
  EM_unif_v     = EM_unif_d / EM_normalise;
  
  all_loglik[0] = sum(log(EM_normalise));
  
  for(int i = 1; i < 100; ++i) {
    
    // M step, update parameters
    EM_exp_sum   = sum(EM_exp);
    EM_norm1_sum = sum(EM_norm1);
    EM_norm2_sum = sum(EM_norm2);
    EM_unif_sum  = sum(EM_unif_v);
    
    EM_parprior[0] = sum( x * EM_exp) / EM_exp_sum;
    EM_parprior[1] = sum( x * EM_norm1) / EM_norm1_sum;
    EM_parprior[2] = sum(pow(x - EM_parprior[1], 2) * EM_norm1) / EM_norm1_sum;
    EM_parprior[3] = sum( x * EM_norm2) / EM_norm2_sum;
    EM_parprior[4] = sum(pow(x - EM_parprior[3], 2) * EM_norm2) / EM_norm2_sum;
    
    // M step, update priors
    EM_parprior[6] = EM_exp_sum   / n;
    EM_parprior[7] = EM_norm1_sum / n;
    EM_parprior[8] = EM_norm2_sum / n;
    EM_parprior[9] = EM_unif_sum  / n;
    
    // M step, restrict parameters
    if (EM_parprior[0] < minexppar) EM_parprior[0] = minexppar;
    if (EM_parprior[6] < 0.05) EM_parprior[6] = 0.05;
    if (EM_parprior[9] > 0.95) EM_parprior[9] = 0.95;
    
    // E step, update posterior
    EM_exp        = (1/EM_parprior[0] * exp(-1/EM_parprior[0] * x)) * EM_parprior[6];
    EM_norm1      = exp(-pow(x - EM_parprior[1], 2) / 2 / EM_parprior[2]) / (2.5066282746 * sqrt(EM_parprior[2])) * EM_parprior[7];
    EM_norm2      = exp(-pow(x - EM_parprior[3], 2) / 2 / EM_parprior[4]) / (2.5066282746 * sqrt(EM_parprior[4])) * EM_parprior[8];
    EM_unif_d     =  1 / (EM_parprior[5]) * EM_parprior[9];
    
    EM_normalise  = EM_exp + EM_norm1 + EM_norm2 + EM_unif_d;
    
    EM_exp        = EM_exp  / EM_normalise;
    EM_norm1      = EM_norm1 / EM_normalise;
    EM_norm2      = EM_norm2 / EM_normalise;
    EM_unif_v     = EM_unif_d / EM_normalise;
    
    all_loglik[i] = sum(log(EM_normalise));
    
    if ((all_loglik[i] - all_loglik[i - 1]) / -all_loglik[i - 1] < tol) break;
    
  }
  
  List final = List::create(Named("pars_priors") = EM_parprior, Named("loglike") = all_loglik);
  
  return  final;    
}   

// [[Rcpp::export]]
List EM_triplemixf(NumericVector x, NumericVector parprior, double tol = 1e-6, double minexppar = 5000) {
  
  int n = x.size();
  NumericVector EM_exp (n);
  NumericVector EM_norm (n);
  NumericVector EM_unif_v (n);
  NumericVector EM_normalise (n);
  NumericVector EM_parprior = parprior;
  NumericVector all_loglik (100);
  double EM_unif_d   = 0;
  double EM_exp_sum  = 0;
  double EM_norm_sum = 0;
  double EM_unif_sum = 0;
  
  // first E step
  
  EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[4];
  EM_norm       = exp(-pow(x - parprior[1], 2) / 2 / parprior[2]) / (2.5066282746 * sqrt(parprior[2])) * parprior[5] * 2;
  EM_unif_d     =  1 / (parprior[3]) * parprior[6];
  
  EM_normalise  = EM_exp + EM_norm + EM_unif_d;
  
  EM_exp        = EM_exp  / EM_normalise;
  EM_norm       = EM_norm / EM_normalise;
  EM_unif_v     = EM_unif_d / EM_normalise;
  
  all_loglik[0] = sum(log(EM_normalise));
  
  EM_parprior[1] = 0;
  
  for(int i = 1; i < 100; ++i) {
    
    // M step, update parameters
    EM_exp_sum  = sum(EM_exp);
    EM_norm_sum = sum(EM_norm);
    EM_unif_sum = sum(EM_unif_v);
    
    EM_parprior[0] = sum( x * EM_exp) / EM_exp_sum;
    EM_parprior[2] = sum(pow(x, 2) * EM_norm) / EM_norm_sum;
    
    // M step, update priors
    EM_parprior[4] = EM_exp_sum  / n;
    EM_parprior[5] = EM_norm_sum / n;
    EM_parprior[6] = EM_unif_sum / n;
    
    // M step, restrict parameters
    if (EM_parprior[0] < minexppar) EM_parprior[0] = minexppar;
    if (EM_parprior[4] < 0.01) EM_parprior[4] = 0.01;
    if (EM_parprior[6] > 0.99) EM_parprior[6] = 0.99;
    
    // E step, update posterior
    EM_exp    = (1/EM_parprior[0] * exp(-1/EM_parprior[0] * x)) * EM_parprior[4];
    EM_norm   = exp(-pow(x - EM_parprior[1], 2) / 2 / EM_parprior[2]) / (2.5066282746 * sqrt(EM_parprior[2])) * EM_parprior[5] * 2;
    EM_unif_d =  1 / (EM_parprior[3]) * EM_parprior[6];
    
    EM_normalise = EM_exp + EM_norm + EM_unif_d;
    
    EM_exp    = EM_exp  / EM_normalise;
    EM_norm   = EM_norm / EM_normalise;
    EM_unif_v = EM_unif_d / EM_normalise;
    
    all_loglik[i] = sum(log(EM_normalise));
    
    if ((all_loglik[i] - all_loglik[i - 1]) / -all_loglik[i - 1] < tol) break;
    
  }
  
  List final = List::create(Named("pars_priors") = EM_parprior, Named("loglike") = all_loglik);
  
  return  final;    
}   

// [[Rcpp::export]]
List EM_quadruplemixf(NumericVector x, NumericVector parprior, double tol = 1e-6, double minexppar = 5000) {
  
  int n = x.size();
  NumericVector EM_exp (n);
  NumericVector EM_norm1 (n);
  NumericVector EM_norm2 (n);
  NumericVector EM_unif_v (n);
  NumericVector EM_normalise (n);
  NumericVector EM_parprior = parprior;
  NumericVector all_loglik (100);
  double EM_unif_d   = 0;
  double EM_exp_sum  = 0;
  double EM_norm1_sum = 0;
  double EM_norm2_sum = 0;
  double EM_unif_sum = 0;
  
  // first E step
  
  EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[6];
  EM_norm1      = exp(-pow(x - parprior[1], 2) / 2 / parprior[2]) / (2.5066282746 * sqrt(parprior[2])) * parprior[7] * 2;
  EM_norm2      = exp(-pow(x - parprior[3], 2) / 2 / parprior[4]) / (2.5066282746 * sqrt(parprior[4])) * parprior[8];
  EM_unif_d     =  1 / (parprior[5]) * parprior[9];
  
  EM_normalise  = EM_exp + EM_norm1 + EM_norm2 + EM_unif_d;
  
  EM_exp        = EM_exp  / EM_normalise;
  EM_norm1      = EM_norm1 / EM_normalise;
  EM_norm2      = EM_norm2 / EM_normalise;
  EM_unif_v     = EM_unif_d / EM_normalise;
  
  all_loglik[0] = sum(log(EM_normalise));
  
  EM_parprior[1] = 0;
  
  for(int i = 1; i < 100; ++i) {
    
    // M step, update parameters
    EM_exp_sum   = sum(EM_exp);
    EM_norm1_sum = sum(EM_norm1);
    EM_norm2_sum = sum(EM_norm2);
    EM_unif_sum  = sum(EM_unif_v);
    
    EM_parprior[0] = sum( x * EM_exp) / EM_exp_sum;
    EM_parprior[2] = sum(pow(x - EM_parprior[1], 2) * EM_norm1) / EM_norm1_sum;
    EM_parprior[3] = sum( x * EM_norm2) / EM_norm2_sum;
    EM_parprior[4] = sum(pow(x - EM_parprior[3], 2) * EM_norm2) / EM_norm2_sum;
    
    // M step, update priors
    EM_parprior[6] = EM_exp_sum   / n;
    EM_parprior[7] = EM_norm1_sum / n;
    EM_parprior[8] = EM_norm2_sum / n;
    EM_parprior[9] = EM_unif_sum  / n;
    
    // M step, restrict parameters
    if (EM_parprior[0] < minexppar) EM_parprior[0] = minexppar;
    if (EM_parprior[6] < 0.01) EM_parprior[6] = 0.01;
    if (EM_parprior[9] > 0.99) EM_parprior[9] = 0.99;
    
    // E step, update posterior
    EM_exp        = (1/EM_parprior[0] * exp(-1/EM_parprior[0] * x)) * EM_parprior[6];
    EM_norm1      = exp(-pow(x - EM_parprior[1], 2) / 2 / EM_parprior[2]) / (2.5066282746 * sqrt(EM_parprior[2])) * EM_parprior[7] * 2;
    EM_norm2      = exp(-pow(x - EM_parprior[3], 2) / 2 / EM_parprior[4]) / (2.5066282746 * sqrt(EM_parprior[4])) * EM_parprior[8];
    EM_unif_d     =  1 / (EM_parprior[5]) * EM_parprior[9];
    
    EM_normalise  = EM_exp + EM_norm1 + EM_norm2 + EM_unif_d;
    
    EM_exp        = EM_exp  / EM_normalise;
    EM_norm1      = EM_norm1 / EM_normalise;
    EM_norm2      = EM_norm2 / EM_normalise;
    EM_unif_v     = EM_unif_d / EM_normalise;
    
    all_loglik[i] = sum(log(EM_normalise));
    
    if ((all_loglik[i] - all_loglik[i - 1]) / -all_loglik[i - 1] < tol) break;
    
  }
  
  List final = List::create(Named("pars_priors") = EM_parprior, Named("loglike") = all_loglik);
  
  return  final;    
}   

// [[Rcpp::export]]
double EM_loglike(NumericVector x, NumericVector parprior) {
  
  int n = x.size();
  NumericVector EM_exp (n);
  NumericVector EM_normalise (n);
  double EM_unif_d   = 0;
  
  if(parprior.size() == 4){
    
    EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[2];
    EM_unif_d     =  1 / (parprior[1]) * parprior[3];
    
    EM_normalise  = EM_exp + EM_unif_d;
  }
  
  if(parprior.size() == 7){
    
    NumericVector EM_norm (n);
    
    EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[4];
    EM_norm       = exp(-pow(x - parprior[1], 2) / 2 / parprior[2]) / (2.5066282746 * sqrt(parprior[2])) * parprior[5];
    EM_unif_d     =  1 / (parprior[3]) * parprior[6];
    
    EM_normalise  = EM_exp + EM_norm + EM_unif_d;
  }
  
  if(parprior.size() == 10){
    
    NumericVector EM_norm1 (n);
    NumericVector EM_norm2 (n);
    
    EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[6];
    EM_norm1      = exp(-pow(x - parprior[1], 2) / 2 / parprior[2]) / (2.5066282746 * sqrt(parprior[2])) * parprior[7];
    EM_norm2      = exp(-pow(x - parprior[3], 2) / 2 / parprior[4]) / (2.5066282746 * sqrt(parprior[4])) * parprior[8];
    EM_unif_d     =  1 / (parprior[5]) * parprior[9];
    
    EM_normalise  = EM_exp + EM_norm1 + EM_norm2 + EM_unif_d;
    
  }
  
  return sum(log(EM_normalise));
}   

// [[Rcpp::export]]
double EM_loglikef(NumericVector x, NumericVector parprior) {
  
  int n = x.size();
  NumericVector EM_exp (n);
  NumericVector EM_normalise (n);
  double EM_unif_d   = 0;
  
  if(parprior.size() == 4){
    
    EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[2];
    EM_unif_d     =  1 / (parprior[1]) * parprior[3];
    
    EM_normalise  = EM_exp + EM_unif_d;
  }
  
  if(parprior.size() == 7){
    
    NumericVector EM_norm (n);
    
    EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[4];
    EM_norm       = exp(-pow(x - parprior[1], 2) / 2 / parprior[2]) / (2.5066282746 * sqrt(parprior[2])) * parprior[5] * 2;
    EM_unif_d     =  1 / (parprior[3]) * parprior[6];
    
    EM_normalise  = EM_exp + EM_norm + EM_unif_d;
  }
  
  if(parprior.size() == 10){
    
    NumericVector EM_norm1 (n);
    NumericVector EM_norm2 (n);
    
    EM_exp        = (1/parprior[0] * exp(-1/parprior[0] * x)) * parprior[6];
    EM_norm1      = exp(-pow(x - parprior[1], 2) / 2 / parprior[2]) / (2.5066282746 * sqrt(parprior[2])) * parprior[7] * 2;
    EM_norm2      = exp(-pow(x - parprior[3], 2) / 2 / parprior[4]) / (2.5066282746 * sqrt(parprior[4])) * parprior[8];
    EM_unif_d     =  1 / (parprior[5]) * parprior[9];
    
    EM_normalise  = EM_exp + EM_norm1 + EM_norm2 + EM_unif_d;
    
  }
  
  return sum(log(EM_normalise));
}   


