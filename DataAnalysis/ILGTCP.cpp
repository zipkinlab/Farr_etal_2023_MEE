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
  DATA_VECTOR( nodes );
  DATA_VECTOR( nobs1 );
  DATA_VECTOR( nobs2 );

  // Projection & weight data
  DATA_VECTOR( weight ); // node weight (dataset 1)
  DATA_SPARSE_MATRIX( A1 ); // integration points
  DATA_SPARSE_MATRIX( A2 ); // observed points (dataset 1)
  DATA_SPARSE_MATRIX( A3 ); // observed points (dataset 2)

  // Intensity design matrix
  DATA_MATRIX( Xmat0 ); // @ nodes
  DATA_MATRIX( Xmat1 ); // @ observations (dataset 1)
  DATA_MATRIX( Xmat2 ); // @ observations (dataset 2)

  // Thinning design matrix
  DATA_MATRIX( Pmat1_1 ); // @ nodes (dataset 1)
  DATA_MATRIX( Pmat2_1 ); // @ observations (dataset 1)

  DATA_MATRIX( Pmat1_2 ); // @ nodes (dataset 1)
  DATA_MATRIX( Pmat2_2 ); // @ observations (dataset 1)

  // SPDE objects
  DATA_STRUCT(spde,spde_t);

  // Fixed effects
  PARAMETER_VECTOR( beta );    // Intensity function intercept & effect parameters
  PARAMETER_VECTOR( alpha_1 );   // Thinning function intercept & effect parameters (dataset1)
  PARAMETER_VECTOR( alpha_2 );   // Thinning function intercept & effect parameters (dataset2)
  PARAMETER( log_kappa );      // Matern parameter

  // Random effects
  PARAMETER_VECTOR( omega ); // Spatial random effect

  Type kappa = exp(log_kappa);

  vector<Type> jnll_comp(5);
  jnll_comp.setZero();

  // Probability of random effects
  SparseMatrix<Type> Q = Q_spde(spde,kappa); // Matern
  jnll_comp(4) = GMRF(Q)( omega );

  // Integration nodes
  vector<Type> nu1 (nodes);
  vector<Type> nu2 (nodes);
  nu1 = weight * exp(Pmat1_1 * alpha_1 + Xmat0 * beta + A1 * omega) / (exp(Pmat1_1 * alpha_1) + 1);
  jnll_comp(3) = nu1.sum();
  nu2 = weight * exp(Pmat1_2 * alpha_2 + Xmat0 * beta + A1 * omega) / (exp(Pmat1_2 * alpha_2) + 1);
  jnll_comp(2) = nu2.sum();

  // Observations (dataset 1)
  vector<Type> mu1 (nobs1);
  mu1 = -1 * (Pmat2_1 * alpha_1 + Xmat1 * beta + A2 * omega - log(exp(Pmat2_1 * alpha_1) + 1));
  jnll_comp(1) = mu1.sum();

  // Observations (dataset 2)
  vector<Type> mu2 (nobs2);
  mu2 = -1 * (Pmat2_2 * alpha_2 + Xmat2 * beta + A3 * omega - log(exp(Pmat2_2 * alpha_2) + 1));
  jnll_comp(0) = mu2.sum();

  // Derived values
  vector<Type> intensity (nodes);
  intensity = exp(Xmat0 * beta + omega);

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( intensity );

  return jnll;
}
