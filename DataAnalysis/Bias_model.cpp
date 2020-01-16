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
  DATA_INTEGER( t_n );
  DATA_INTEGER( p_n );
  DATA_IVECTOR( tp_i );
  DATA_IVECTOR( tp_k );


  // Projection & weight data
  DATA_VECTOR( weight ); // node weight
  DATA_SPARSE_MATRIX( A1 ); // projection matrix spring1
  DATA_SPARSE_MATRIX( A2 ); // projection matrix spring2
  DATA_SPARSE_MATRIX( A3 ); // projection matrix summer

  // SPDE objects
  DATA_STRUCT(spde_sp,spde_t);
  DATA_STRUCT(spde_su,spde_t);

  // Fixed effects
  PARAMETER( beta0 );   // Intensity function intercept by year
  PARAMETER( delta1 );
  PARAMETER( delta2 );
  PARAMETER( gamma1 );
  PARAMETER( gamma2 );
  PARAMETER( log_kappa_sp1 );      // Scale parameter spring 1
  PARAMETER( log_kappa_sp2 );      // Scale parameter spring 2
  PARAMETER( log_kappa_su );      // Scale parameter summer
  PARAMETER( log_tau_sp1 );      // Precision parameter for spatial effect spring 1
  PARAMETER( log_tau_sp2 );      // Precision parameter for spatial effect spring 2
  PARAMETER( log_tau_su );      // Precision parameter for spatial effect summer

  // Random effects
  PARAMETER_VECTOR( omega_sp1 ); // Spatial random effect spring 1
  PARAMETER_VECTOR( omega_sp2 ); // Spatial random effect spring 2
  PARAMETER_VECTOR( omega_su ); // Spatial random effect summer

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

  vector<Type> jnll_comp(9);
  jnll_comp.setZero();

  // Probability of random effects
  SparseMatrix<Type> Qsp1 = Q_spde(spde_sp,kappa_sp1);
  SparseMatrix<Type> Qsp2 = Q_spde(spde_sp,kappa_sp2);
  SparseMatrix<Type> Qsu = Q_spde(spde_su,kappa_su);
  jnll_comp(0) += GMRF(Qsp1)( omega_sp1 );
  jnll_comp(1) += GMRF(Qsp2)( omega_sp2 );
  jnll_comp(2) += GMRF(Qsu)( omega_su );

  // Holding values
  vector<Type> Omega_sp1(nodes_sp);
  vector<Type> Omega_sp2(nodes_sp);
  vector<Type> Omega_su(nodes_su);
  //matrix<Type> Omega(nodes, ??); // This should be 2
  //matrix<Type> Epsilon(nodes, t_n);
  vector<Type> omg_sp1(nobs_sp1);
  vector<Type> omg_sp2(nobs_sp2);
  vector<Type> omg_su(nobs_su);

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
      jnll_comp(3) += weight(j) * exp(beta(tp_k(j)) + Omega_sp1(k));
      jnll_comp(4) += weight(i) * exp(beta(tp_k(i)) + Omega_sp2(k));
    }

    for(int k=0; k<nodes_su; k++, g++){
      jnll_comp(5) += weight(g) * exp(beta(tp_k(g)) + Omega_su(k));
    }
  }

  // Intensity function @ observations
  int h=0;
  for(int i=0; i<nobs_sp1; i++, h++){
    jnll_comp(6) -= beta(tp_i(h)) + omg_sp1(i);
  }

  for(int i=0; i<nobs_sp2; i++, h++){
    jnll_comp(7) -= beta(tp_i(h)) + omg_sp2(i);
  }

  for(int i=0; i<nobs_su; i++, h++){
    jnll_comp(8) -= beta(tp_i(h)) + omg_su(i);
  }

  // Joint NLL
  Type jnll = jnll_comp.sum();

  // Reporting
  REPORT( jnll_comp );
  REPORT( jnll );

  return jnll;
}
