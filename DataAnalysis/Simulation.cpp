#include <TMB.hpp>
// Integrated Log Gaussian Cox Process
template <class Type>
Type objective_function<Type>::operator() ()
{
  
  // objective function -- joint negative log-likelihood
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // DATA //
  
  // Indices
  DATA_INTEGER( nodes );  // Number of nodes 
  DATA_INTEGER( nobs );   // Number of presence-only observations
  DATA_INTEGER( ncount ); // Number of structured sampling sites
  DATA_INTEGER( t_n );    // Number of stages
  DATA_IVECTOR( t_i );    // Stage ID for data
  
  // Projection & weight data
  DATA_VECTOR( weight );   // Node weight
  DATA_VECTOR( area );     // Area swept for counts
  DATA_SPARSE_MATRIX( A ); // Projection matrix
  
  // Count dataset
  DATA_VECTOR( counts ); // Number of counts per site
  
  // Covaraite dataset
  DATA_VECTOR( nBias ); //Sampling bias @ nodes
  DATA_MATRIX( nCov );  //Environmental covariate @ nodes
  
  DATA_VECTOR( Bias ); //Sampling bias @ observations
  DATA_VECTOR( Cov ); //Environmental covariate @ data
  
  // SPDE objects
  DATA_STRUCT(spde,spde_t); // Sparse matrix for Matern covariance structure
  
  // Prediction
  DATA_INTEGER( npred );        // Number of pixels for prediction
  DATA_MATRIX( predX );         // Environmental covariate for prediction
  DATA_SPARSE_MATRIX( Apred );  // Projection matrix for prediction
  
  // PARAMETERS //
  
  // Fixed effects
  PARAMETER( beta0 );     // Baseline population density
  PARAMETER( beta1 );     // Effect of environmental covariate
  PARAMETER( delta1 );    // Effect of stage (population change)
  PARAMETER( alpha0 );    // Thinning function intercept
  PARAMETER( alpha1 );    // Effect of sampling bias
  PARAMETER( log_kappa ); // Scale parameter of Matern covariance
  PARAMETER( log_tau_O ); // Precision parameter of Matern covariance
  
  // Random effects
  PARAMETER_VECTOR( omega ); // Spatial random effect
  
  // Population density at each stage
  vector<Type> beta(t_n);
  beta(0) = beta0;
  beta(1) = beta0 + delta1;

  // Derived parameters
  Type kappa = exp(log_kappa);
  Type tau_O = exp(log_tau_O);
  Type range = sqrt(8)/kappa;
  Type sigma_O = 1/sqrt(4*PI*tau_O*tau_O*kappa*kappa);
  
  vector<Type> jnll_comp(4);
  jnll_comp.setZero();
  
  // Probability of random effects
  SparseMatrix<Type> Q = Q_spde(spde,kappa); // Matern covariance (see R_inla namespace)
  jnll_comp(0) += GMRF(Q)( omega );
  
  // Holding values
  vector<Type> Omega(nodes);
  vector<Type> omg(nobs + ncount);
  
  // Transform GMRFs
  for(int k=0; k<nodes; k++){
    Omega(k) = omega(k) / tau_O;
  }
  
  // Project GMRFs
  omg = A * Omega; // Project omega to points
  
  // Intensity function @ nodes
  for(int t=0; t<t_n; t++){
    for(int k=0; k<nodes; k++){
      jnll_comp(1) += weight(k) * exp(alpha0 + alpha1 * nBias(k) + beta(t) + beta1 * nCov(k,t)+ Omega(k)) / (exp(alpha0 + alpha1 * nBias(k)) + 1); // Integration nodes
    }
  }
  
  // Intensity function @ presence-only data
  for(int i=0; i<nobs; i++){
    jnll_comp(2) -= alpha0 + alpha1 * Bias(i) + beta(t_i(i)) + beta1 * Cov(i) + omg(i) - log(exp(alpha0 + alpha1 * Bias(i)) + 1); // Observation points
  }
  
  // Intensity function @ counts
  vector<Type> lambda(ncount);
  for(int j=0; j<ncount; j++){
    lambda(j) = area(j) * exp(beta0 + beta1 * Cov(j+nobs) + omg(j+nobs));
    jnll_comp(3) -= dpois(counts(j), lambda(j), true);
  }
  
  // Prediction
  vector<Type> predO(npred);
  predO = Apred * Omega;
  
  vector<Type> pred(t_n);
  vector<Type> Npred(t_n);
  
  for(int t=0; t<t_n; t++){
    for(int g=0; g<npred; g++){
      pred(t) += exp(beta(t) + beta1 * predX(g,t) + predO(g));
    }
    Npred(t) = pred(t)/npred * 9;
  }
  
  
  // Joint NLL
  Type jnll = jnll_comp.sum();
  
  // Reporting
  REPORT( kappa );
  REPORT( tau_O );
  REPORT( sigma_O );
  REPORT( range );
  REPORT( Npred );
  ADREPORT( Npred );
  ADREPORT( kappa );
  ADREPORT( tau_O );
  ADREPORT( sigma_O );
  
  REPORT( jnll_comp );
  REPORT( jnll );
  
  return jnll;
}