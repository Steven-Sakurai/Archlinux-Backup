/*
	Calculate the initial filter moments $M_{k, 0}^{(l)}$
	In this file, we set 
		k = [0, 20] (22 nodes)
		l = [0, 3] (cubic piecewise poly)
*/


#ifndef _INITIALMOMENT
#define _INITIALMOMENT

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;
typedef vector< double > stdVec;
typedef vector< vector<double> > stdMat;

#include <boost/math/special_functions/gamma.hpp>
using boost::math::gamma_p;

#include "Parameter.h" 

// the stdMat should be sth like `stdMat initialMoments(4, stdVec(21, 0));`
void initialMoments(Parameter par, stdMat initialMoments, bool output) {
    double mu = par.mu;
    double kappa = par.kappa;
    double theta = par.theta;
    double xi = par.xi;
    double rho = par.rho;
	// scale parameter of gamma dis
	double beta = xi*xi/(2*kappa);
	// taking log of the mean and the standard deviation of gamma dis
	double stLogMean = log(theta);
	double stLogSD = log(sqrt(theta*xi*xi/(2*kappa)));
	// create the nodes we want to interpolate
	stdVec interNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35,
		 -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};
	for(int i = 0; i < 22; ++i) 
		interNodes[i] = exp(stLogMean - interNodes[i]*stLogSD) / (xi*xi/(2*kappa));
	// shape parameter of gamma dis
	double shape = 2*kappa*theta/(xi*xi);
	// calculate using the beautiful property of gamma distribution	
	for(int i = 0; i < 21; ++i)
		initialMoments[0][i] = gamma_p(shape, interNodes[i+1]) - gamma_p(shape, interNodes[i]);
	for(int i = 0; i < 21; ++i)
		initialMoments[1][i] = shape*beta*(gamma_p(shape+1, interNodes[i+1]) - gamma_p(shape+1, interNodes[i]));
	for(int i = 0; i < 21; ++i)
		initialMoments[2][i] = shape*(shape+1)*beta*beta * (gamma_p(shape+2, interNodes[i+1]) - gamma_p(shape+2, interNodes[i]));
	for(int i = 0; i < 21; ++i)
		initialMoments[3][i] = shape*(shape+1)*(shape+2)*beta*beta*beta * (gamma_p(shape+3, interNodes[i+1]) - gamma_p(shape+3, interNodes[i]));
	// end
	// output part
	if(output == true) {
		cout << "0th order basis:" << setprecision(16) << endl;
		for(int i = 0; i < 21; ++i) {
			cout << initialMoments[0][i] << '\t';
			if((i+1)%5 == 0)
				cout << endl;
		}
		cout << endl;
		cout << "1th order basis:" << endl;
		for(int i = 0; i < 21; ++i) {
			cout << initialMoments[1][i] << '\t';
			if((i+1)%5 == 0)
				cout << endl;
		}
		cout << endl;
		cout << "2th order basis:" << endl;
		for(int i = 0; i < 21; ++i) {
			cout << initialMoments[2][i] << '\t';
			if((i+1)%5 == 0)
				cout << endl;
		}
		cout << endl;
		cout << "3th order basis:" << endl;
		for(int i = 0; i < 21; ++i) {
			cout << initialMoments[3][i] << '\t';
			if((i+1)%5 == 0)
				cout << endl;
		}
	}
}




#endif