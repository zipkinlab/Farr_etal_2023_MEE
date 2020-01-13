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
  DATA_INTEGER( nodes_sp );
  DATA_INTEGER( nodes_su );
  DATA_INTEGER( nobs_sp1 );
  DATA_INTEGER( nobs_sp2 );
  DATA_INTEGER( nobs_su );
  DATA_INTEGER( ncount );
  DATA_INTEGER( ncount_sp1 );
  DATA_INTEGER( ncount_sp2 );
  DATA_INTEGER( ncount_su );
  DATA_INTEGER( t_n );
  DATA_INTEGER( p_n );
  DATA_IVECTOR( tp_i );
  DATA_IVECTOR( tp_k );
  // Projection & weight data
  DATA_VECTOR( weight ); // node weight
  DATA_VECTOR( area ); // area sweept for counts
  DATA_SPARSE_MATRIX( A1 ); // projection matrix spring1
  DATA_SPARSE_MATRIX( A2 ); // projection matrix spring2
  DATA_SPARSE_MATRIX( A3 ); // projection matrix summer

  //Count dataset
  DATA_VECTOR( counts );

  // Covaraite dataset
  DATA_VECTOR( nBias ); //Bias @ nodes
  DATA_VECTOR( nCov ); //Covariate @ nodes

  DATA_VECTOR( Bias ); //Bias @ observations
  DATA_VECTOR( Cov ); //Covariate @ observations

  // SPDE objects
  DATA_STRUCT(spde_sp,spde_t);
  DATA_STRUCT(spde_su,spde_t);

  // Fixed effects
  PARAMETER( beta0 );   // Intensity function intercept by year
  PARAMETER( beta1 );          // Effect of NDVI on intensity
  PARAMETER( delta1 );
  PARAMETER( delta2 );
  PARAMETER( gamma1 );
  PARAMETER( gamma2 );
  PARAMETER( alpha0 );         // Thinning function intercept
  PARAMETER( alpha1 );         // Effect of sampling bias
  PARAMETER( log_kappa_sp1 );      // Scale parameter spring 1
  PARAMETER( log_kappa_sp2 );      // Scale parameter spring 2
  PARAMETER( log_kappa_su );      // Scale parameter summer
  PARAMETER( log_tau_sp1 );      // Precision parameter for spatial effect spring 1
  PARAMETER( log_tau_sp2 );      // Precision parameter for spatial effect spring 2
  PARAMETER( log_tau_su );      // Precision parameter for spatial effect summer
  //PARAMETER( log_tau_E );      // Precision parameter for spatial effect by time

  // Random effects
  PARAMETER_VECTOR( omega_sp1 ); // Spatial random effect spring 1
  PARAMETER_VECTOR( omega_sp2 ); // Spatial random effect spring 2
  PARAMETER_VECTOR( omega_su ); // Spatial random effect summer
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

  Type kappa_sp1 = exp(log_kappa_sp1);
  Type kappa_sp2 = exp(log_kappa_sp2);
  Type kappa_su = exp(log_kappa_su);
  Type tau_sp1 = exp(log_tau_sp1);
  Type tau_sp2 = exp(log_tau_sp2);
  Type tau_su = exp(log_tau_su);
  //Type tau_E = exp(log_tau_E);
  //Type range = sqrt(8)/kappa;
  //Type sigma_O = 1/sqrt(4*PI*tau_O*tau_O*kappa*kappa);
  //Type sigma_E = 1/sqrt(4*PI*tau_E*tau_E*kappa*kappa);

  //vector<Type> jnll_comp(5);
  vector<Type> jnll_comp(12);
  jnll_comp.setZero();

  // Probability of random effects
  SparseMatrix<Type> Qsp1 = Q_spde(spde_sp,kappa_sp1);
  SparseMatrix<Type> Qsp2 = Q_spde(spde_sp,kappa_sp2);
  SparseMatrix<Type> Qsu = Q_spde(spde_su,kappa_su);
  jnll_comp(0) += GMRF(Qsp1)( omega_sp1 );
  jnll_comp(1) += GMRF(Qsp2)( omega_sp2 );
  jnll_comp(2) += GMRF(Qsu)( omega_su );
  //for(int t=0; t<t_n; t++){
  //  jnll_comp(1) += GMRF(Q)( epsilon.col(t) );
  //}

  // Holding values
  vector<Type> Omega_sp1(nodes_sp);
  vector<Type> Omega_sp2(nodes_sp);
  vector<Type> Omega_su(nodes_su);
  //matrix<Type> Omega(nodes, ??); // This should be 2
  //matrix<Type> Epsilon(nodes, t_n);
  vector<Type> omg_sp1(nobs_sp1 + ncount_sp1);
  vector<Type> omg_sp2(nobs_sp2 + ncount_sp2);
  vector<Type> omg_su(nobs_su + ncount_su);
  //vector<Type> omg(nobs + ncount);
  //matrix<Type> eps((nobs + ncount), t_n);

  // Transform GMRFs

  //Spring
  for(int k=0; k<nodes_sp; k++){
    Omega_sp1(k) = omega_sp1(k) / tau_sp1;
    Omega_sp2(k) = omega_sp2(k) / tau_sp2;
  }
  //Summer
  for(int k=0; k<nodes_su; k++){
    Omega_su(k) = omega_su(k) / tau_su;
  }

  // Project GMRFs
  omg_sp1 = A1 * Omega_sp1; //Spring 1; length 257
  omg_sp2 = A2 * Omega_sp2; //Spring 2; length 303
  omg_su = A3 * Omega_su; //Summer; length 825

  // Intensity function @ nodes
  for(int t=0, j=0, i=nodes_sp*t_n, g=nodes_sp*t_n*2; t<t_n; t++){
    for(int k=0; k<nodes_sp; k++, j++, i++){
      jnll_comp(3) += weight(j) * exp(alpha0 + alpha1 * nBias(j) + beta(tp_k(j)) + beta1 * nCov(j) + Omega_sp1(k)) / (exp(alpha0 + alpha1 * nBias(j)) + 1);
      jnll_comp(4) += weight(i) * exp(alpha0 + alpha1 * nBias(i) + beta(tp_k(i)) + beta1 * nCov(i) + Omega_sp2(k)) / (exp(alpha0 + alpha1 * nBias(i)) + 1);
    }

    for(int k=0; k<nodes_su; k++, g++){
      jnll_comp(5) += weight(g) * exp(alpha0 + alpha1 * nBias(g) + beta(tp_k(g)) + beta1 * nCov(g) + Omega_su(k)) / (exp(alpha0 + alpha1 * nBias(g)) + 1);
    }
  }
  // Intensity function @ points
  int j=0;
  for(int i=0; i<nobs_sp1; i++, j++){
    jnll_comp(6) -= alpha0 + alpha1 * Bias(j) + beta(tp_i(j)) + beta1 * Cov(j) + omg_sp1(i) - log(exp(alpha0 + alpha1 * Bias(j)) + 1);
  }

  for(int i=0; i<nobs_sp2; i++, j++){
    jnll_comp(7) -= alpha0 + alpha1 * Bias(j) + beta(tp_i(j)) + beta1 * Cov(j) + omg_sp2(i) - log(exp(alpha0 + alpha1 * Bias(j)) + 1);
  }

  for(int i=0; i<nobs_su; i++, j++){
    jnll_comp(8) -= alpha0 + alpha1 * Bias(j) + beta(tp_i(j)) + beta1 * Cov(j) + omg_su(i) - log(exp(alpha0 + alpha1 * Bias(j)) + 1);
  }

  // Intensity function @ counts
  vector<Type> lambda(ncount);
  int g=0;
  for(int i=0; i<ncount_sp1; i++, j++, g++){
    lambda(g) = area(g) * exp(beta(tp_i(j)) + beta1 * Cov(j) + omg_sp1(i+nobs_sp1));
    jnll_comp(9) -= dpois(counts(g), lambda(g), true);
  }

  for(int i=0; i<ncount_sp2; i++, j++, g++){
    lambda(g) = area(g) * exp(beta(tp_i(j)) + beta1 * Cov(j) + omg_sp2(i+nobs_sp2));
    jnll_comp(10) -= dpois(counts(g), lambda(g), true);
  }

  for(int i=0; i<ncount_su; i++, j++, g++){
    lambda(g) = area(g) * exp(beta(tp_i(j)) + beta1 * Cov(j) + omg_su(i+nobs_su));
    jnll_comp(11) -= dpois(counts(g), lambda(g), true);
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
