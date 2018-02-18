// first we need a class to store parameter: mu, kappa, theta, xi, rho
// call the object `class Parameter par`

// N stands for the number of obs of x = log(S)

// Use 21 basis function at each stage, i.e. 21 intervals, 22 nodes for log(S)
// every basis is cubic, need 4 parameters



/*
	need to get Yobs, i.e. simulation
*/

/*
	need to implement BasisCoefficients function
	to get \alpha and \beta
	I think we should use armadillo in this case
*/


#include "InitialMoments.hpp"



double loglike(Parameter Par, stdVec YObs, int N, int n) {
	mu = Par.mu;
	kappa = Par.kappa;
	theta = Par.theta;
	xi = Par.xi;
	rho = Par.rho;
	
}













}