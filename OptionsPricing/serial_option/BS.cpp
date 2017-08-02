/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 * 
 * This file contains routines to parallely compute the call and 
 * put price of an European option and Asian option.
 * 
 * Ruby Zhang
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 

#include <algorithm>    // Needed for the "max" function
#include <cmath>
#include <iostream>
#include <omp.h>


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 A simple implementation of the Box-Muller algorithm, used to 
generate gaussian random numbers; necessary for the Monte Carlo 
method below. */

double gaussian_box_muller(unsigned int& seed) {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  // Continue generating two uniform random variables
  // until the square of their "euclidean distance" 
  // is less than unity
  do {
    x = 2.0 * rand_r(&seed) / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand_r(&seed) / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Pricing a European vanilla call option with a Monte Carlo method

double monte_carlo_call_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T) {
  double S_adjust = S * exp(T*(r-0.5*v*v));
  double S_cur = 0.0;
  double payoff_sum = 0.0;

  #pragma omp parallel for reduction(+:payoff_sum)
  for (unsigned int i=0; i<num_sims; i++) {
    unsigned int seed = i*omp_get_thread_num();
    double gauss_bm = gaussian_box_muller(seed);
    S_cur = S_adjust * exp(sqrt(v*v*T)*gauss_bm);
    payoff_sum += std::max(S_cur - K, 0.0);
  }

  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Pricing a European vanilla put option with a Monte Carlo method

double monte_carlo_put_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T) {
  double S_adjust = S * exp(T*(r-0.5*v*v));
  double S_cur = 0.0;
  double payoff_sum = 0.0;

  #pragma omp parallel for reduction(+:payoff_sum)
  for (unsigned int i=0; i<num_sims; i++) {
    unsigned int seed = i*omp_get_thread_num();
    double gauss_bm = gaussian_box_muller(seed);
    S_cur = S_adjust * exp(sqrt(v*v*T)*gauss_bm);
    payoff_sum += std::max(K - S_cur, 0.0);
  }

  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Pricing an Asian call option with a Monte Carlo method

double asian_monte_carlo_call_price(const int& num_sims, const int& periods, const double& S, const double& K, const double& r, const double& v, const double& T) {
  double S_prev = S;
  double S_cur = 0.0;
  double payoff_sum = 0.0;
  double S_bar = 0.0;
  double interval = T/periods;
  double adjust = exp((r-0.5*v*v)*interval);
  
  #pragma omp parallel for firstprivate(S_bar,S_prev) reduction(+:payoff_sum)
  for (unsigned int i=0; i<num_sims; i++) {
    S_bar = 0.0;
    S_prev = S;
    for (unsigned int j=0; j<periods; j++){
    	unsigned int seed = (i*periods+j)*omp_get_thread_num();
    	double gauss_bm = gaussian_box_muller(seed);
    	S_cur = S_prev*adjust*exp(sqrt(v*v*interval)*gauss_bm);
	S_prev = S_cur;
	S_bar += S_cur;	
    }
    S_bar = S_bar/static_cast<double>(periods);
    payoff_sum += std::max(S_bar - K, 0.0);
  }
  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Pricing an Asian put option with a Monte Carlo method

double asian_monte_carlo_put_price(const int& num_sims, const int& periods, const double& S, const double& K, const double& r, const double& v, const double& T) {
  double S_prev = S;
  double S_cur = 0.0;
  double payoff_sum = 0.0;
  double S_bar = 0.0;
  double interval = T/periods;
  double adjust = exp((r-0.5*v*v)*interval);

  #pragma omp parallel for firstprivate(S_bar,S_prev) reduction(+:payoff_sum) 
  for (unsigned int i=0; i<num_sims; i++) {
    S_bar = 0.0;
    S_prev = S;
    for (unsigned int j=0; j<periods; j++){
    	unsigned int seed = (i*periods+j)*omp_get_thread_num();
    	double gauss_bm = gaussian_box_muller(seed);
    	S_cur = S_prev*adjust*exp(sqrt(v*v*interval)*gauss_bm);
	S_prev = S_cur;
	S_bar += S_cur;
    }
    S_bar = S_bar/static_cast<double>(periods);
    payoff_sum += std::max(K - S_bar, 0.0);
  }

  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char **argv) {

  // Parameters                                                                             
  int num_sims = 1000000;   // Number of simulated asset paths                                                       
  double S = 100.0;  // Option price                                                                                  
  double K = 100.0;  // Strike price                                                                                  
  double r = 0.05;   // Risk-free rate (5%)                                                                           
  double v = 0.2;    // Volatility of the underlying (20%)                                                            
  double T = 1.0;    // One year until expiry

  double periods = 12; // Periods of Asian option

  // Then we calculate the call/put values via Monte Carlo                                                                          
  double call1 = monte_carlo_call_price(num_sims, S, K, r, v, T);
  double put1 = monte_carlo_put_price(num_sims, S, K, r, v, T);

  double call2 = asian_monte_carlo_call_price(num_sims, periods, S, K, r, v, T);
  double put2 = asian_monte_carlo_put_price(num_sims, periods, S, K, r, v, T);

  // Finally we output the parameters and prices                                                                      
  std::cout << "Number of Paths: " << num_sims << std::endl;
  std::cout << "Number of Periods: " <<  periods << std::endl;
  std::cout << "Underlying:      " << S << std::endl;
  std::cout << "Strike:          " << K << std::endl;
  std::cout << "Risk-Free Rate:  " << r << std::endl;
  std::cout << "Volatility:      " << v << std::endl;
  std::cout << "Maturity:        " << T << std::endl;

  std::cout << "European Call Price:   " << call1 << std::endl;
  std::cout << "European Put Price:    " << put1 << std::endl;

  std::cout << "Asian Call Price:   " << call2 << std::endl;
  std::cout << "Asian Put Price:    " << put2 << std::endl;

  return 0;
}
