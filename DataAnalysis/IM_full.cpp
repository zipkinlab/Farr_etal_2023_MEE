#include <TMB.hpp>
// Integrated Log Gaussian Cox Process
template <class Type>
Type objective_function<Type>::operator() ()
{

  // objective function -- joint negative log-likelihood
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  //Indices
  DATA_INTEGER( nodes );
  DATA_INTEGER( nobs );
  DATA_INTEGER( ncount );
  DATA_INTEGER( t_n );
  DATA_INTEGER( p_n );
  DATA_IVECTOR( tp_i );
  DATA_IVECTOR( tp_k );
  // Projection & weight data
  DATA_VECTOR( weight ); // node weight
  DATA_VECTOR( area ); // area sweept for counts
  //DATA_SPARSE_MATRIX( A1 ); // projection matrix spring
  //DATA_SPARSE_MATRIX( A2 ); // projection matrix summer

  //Count dataset
  DATA_VECTOR( counts );

  // Covaraite dataset
  DATA_VECTOR( nBias ); //Bias @ nodes
  DATA_VECTOR( nCov ); //Covariate @ nodes

  DATA_VECTOR( Bias ); //Bias @ observations
  DATA_VECTOR( Cov ); //Covariate @ observations

  // SPDE objects
  //DATA_STRUCT(spde,spde_t);

  // Fixed effects
  PARAMETER( beta0 );   // Intensity function intercept by year
  PARAMETER( beta1 );          // Effect of NDVI on intensity
  PARAMETER( delta1 );
  PARAMETER( delta2 );
  PARAMETER( gamma1 );
  PARAMETER( gamma2 );
  PARAMETER( alpha0 );         // Thinning function intercept
  PARAMETER( alpha1 );         // Effect of sampling bias
  //PARAMETER( log_kappa );      // Scale parameter
  //PARAMETER( log_tau_O );      // Precision parameter for spatial effect
  //PARAMETER( log_tau_E );      // Precision parameter for spatial effect by time

  // Random effects
  //PARAMETER_VECTOR( omega ); // Spatial random effect
  //PARAMETER_ARRAY( epsilon ); // Spatial process variatin

  vector<Type> beta(t_n * p_n);
  beta(0) = beta0 + delta1 + delta2;
  beta(1) = beta0 + delta1;
  beta(2) = beta0;
  beta(3) = beta0 + gamma1 + delta1 + delta2;
  beta(4) = beta0 + gamma1 + delta1;
  beta(5) = beta0 + gamma1;
  beta(6) = beta0 + gamma1 + gamma2 + delta1 + delta2;
  beta(7) = beta0 + gamma1 + gamma2 + delta1;
  beta(8) = beta0 + gamma1 + gamma2;

  //Type kappa = exp(log_kappa);
  //Type tau_O = exp(log_tau_O);
  //Type tau_E = exp(log_tau_E);
  //Type range = sqrt(8)/kappa;
  //Type sigma_O = 1/sqrt(4*PI*tau_O*tau_O*kappa*kappa);
  //Type sigma_E = 1/sqrt(4*PI*tau_E*tau_E*kappa*kappa);

  //vector<Type> jnll_comp(5);
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();

  // Probability of random effects
  //SparseMatrix<Type> Q = Q_spde(spde,kappa); // Matern covariance (see R_inla namespace)
  //jnll_comp(0) += GMRF(Q)( omega );
  //for(int t=0; t<t_n; t++){
  //  jnll_comp(1) += GMRF(Q)( epsilon.col(t) );
  //}

  // Holding values
  //vector<Type> Omega(nodes);
  //matrix<Type> Omega(nodes, ??); // This should be 2
  //matrix<Type> Epsilon(nodes, t_n);
  //vector<Type> omg(nobs + ncount);
  //matrix<Type> eps((nobs + ncount), t_n);

  // Transform GMRFs

  //Spring
  //for(int k=0; k<nodes_sp; k++){
  //  Omega_sp(k) = omega_sp(k) / tau_sp;
  //}
  //Summer
  //for(int k=0; k<nodes_su; k++){
  //  Omega_su(k) = omega_su(k) / tau_su;
  //}

  //Omega(k, 1) = omega

  // Project GMRFs
  //omg_sp = A1 * Omega_sp; //Spring
  //omg_su = A2 * Omega_su; //Summer

  // Intensity function @ nodes
  for(int k=0; k<nodes; k++){
    //jnll_comp(2) += weight(k) * exp(alpha0 + alpha1 * nBias(k) + beta(t_k(k)) + beta1 * nCov(k) + Omega(k, t_k(k))) / (exp(alpha0 + alpha1 * nBias(k)) + 1);
    jnll_comp(0) += weight(k) * exp(alpha0 + alpha1 * nBias(k) + beta(tp_k(k)) + beta1 * nCov(k)) / (exp(alpha0 + alpha1 * nBias(k)) + 1);
  }

  // Intensity function @ points
  for(int i=0; i<nobs; i++){
    //jnll_comp(3) -= alpha0 + alpha1 * Bias(i) + beta(t_i(i)) + beta1 * Cov(i) + omg(i) + eps(i, t_i(i)) - log(exp(alpha0 + alpha1 * Bias(i)) + 1); // Observation points
    jnll_comp(1) -= alpha0 + alpha1 * Bias(i) + beta(tp_i(i)) + beta1 * Cov(i) - log(exp(alpha0 + alpha1 * Bias(i)) + 1); // Observation points
  }

  // Intensity function @ counts
  vector<Type> lambda(ncount);
  for(int j=0; j<ncount; j++){
    //lambda(j) = area(j) * exp(beta0 + beta1 * Cov(j+nobs) + omg(j+nobs) + eps((j+nobs), 1));
    lambda(j) = area(j) * exp(beta(tp_i(j+nobs)) + beta1 * Cov(j+nobs));
    jnll_comp(2) -= dpois(counts(j), lambda(j), true);
  }

  // Joint NLL
  Type jnll = jnll_comp.sum();

  // Reporting
  //REPORT( tau_O );
  //REPORT( sigma_O );
  //REPORT( tau_E );
  //REPORT( sigma_E );
  //REPORT( range );

  REPORT( jnll_comp );
  REPORT( jnll );

  return jnll;
}
