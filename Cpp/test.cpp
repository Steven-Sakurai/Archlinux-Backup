#include <boost/math/special_functions/gamma.hpp>
using boost::math::gamma_p;

#include <omp.h>

#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <iomanip>

using namespace std;
typedef vector< double > stdVec;
typedef vector< vector< double > > stdMat;

#define Pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

#define Delta 0.003968254 // 1/252
#define nbasis 21
const stdVec myNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};

#define njumps 1
/*
	Transition density parameter
	0 stands for no jump
	1 stands for 1 jump
*/

inline double tran_y_mean_0(const stdVec& par, double y0, double v0) {
	return y0 + (par[0]-0.5*v0) * Delta;
}

inline double tran_y_var_0(const stdVec& par, double v0) {
	return v0*Delta;
}

inline double tran_y_mean_n(const stdVec& par, double y0, double v0, int n) {
	return y0 + (par[0]-0.5*v0)*Delta + par[6]*n;
}

inline double tran_y_var_n(const stdVec& par, double v0, int n) {
	return v0*Delta + par[7]*n;
}


inline double tran_v_mean(const stdVec& par, double v0) {
	return (v0 + par[1]*(par[2] - v0)*Delta);
}

inline double tran_v_var(const stdVec& par, double v0) {
	return pow(par[3], 2)*v0*Delta;
}

/*
	transition density decomposition: p(y, v) = p(y) * p(v|y) 
*/

inline double cond_v_var(const stdVec& par, double v0)
{
	return tran_v_var(par, v0) * (1.0 - pow(par[4], 2));
}

inline double cond_v_mean_0(const stdVec& par, double y, double y0, double v0)
{
	return tran_v_mean(par, v0) + sqrt(tran_v_var(par, v0)/tran_y_var_0(par, v0))*par[4]*(y - tran_y_mean_0(par, y0, v0));
}

inline double cond_v_mean_n(const stdVec& par, double y, double y0, double v0, int n)
{
	return tran_v_mean(par, v0) + sqrt(tran_v_var(par, v0)/tran_y_var_n(par, v0, n))*par[4]*(y - tran_y_mean_n(par, y0, v0, n));
}

inline double prePy_0(const stdVec& par, double y, double y0, double v0) {
	double mu0 = tran_y_mean_0(par, y0, v0);
	double var0 = tran_y_var_0(par, v0);
	return exp(-(y-mu0)*(y-mu0) / (2.0*var0)) / sqrt(2*Pi*var0);
}

inline double prePy_n(const stdVec& par, double y, double y0, double v0, int n) {
	double mu1 = tran_y_mean_n(par, y0, v0, n);
	double var1 = tran_y_var_n(par, v0, n);
	return exp(-(y- mu1)*(y-mu1)/(2.0*var1)) / sqrt(2*Pi*var1);
}

// p(v|y) conditional distribution
inline double pdf_c_0(const stdVec& par, double y, double y0, double v, double v0) {
	double mu = cond_v_mean_0(par, y, y0, v0);
	double var = cond_v_var(par, v0);
	return ( exp(-(v-mu)*(v-mu) / (2.0*var) ) / sqrt(2*Pi*var) );
}

inline double pdf_c_n(const stdVec& par, double y, double y0, double v, double v0, int n) {
	double mu = cond_v_mean_n(par, y, y0, v0, n);
	double var = cond_v_var(par, v0);
	return ( exp(-(v-mu)*(v-mu) / (2.0*var) ) / sqrt(2*Pi*var) );
}

inline double cdf_c_0(const stdVec& par, double y, double y0, double v, double v0) {
	double mu = cond_v_mean_0(par, y, y0, v0);
	double var = cond_v_var(par, v0);
	return 1.0 - 0.5 * erfc((v-mu)/sqrt(2.0*var));
}

inline double cdf_c_n(const stdVec& par, double y, double y0, double v, double v0, int n) {
	double mu = cond_v_mean_n(par, y, y0, v0, n);
	double var = cond_v_var(par, v0);
	return 1.0 - 0.5 * erfc((v-mu)/sqrt(2.0*var));
}


// the four integral B_k
// 0 for no jump, 1 for 1 jump
inline double Integ0_0(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	return prePy_0(par, y, y0, v0) * (cdf_c_0(par, y, y0, v2, v0) - cdf_c_0(par, y, y0, v1, v0));
}

inline double Integ1_0(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	return Integ0_0(par, y, y0, v1, v2, v0)*cond_v_mean_0(par, y, y0, v0)  - 
		cond_v_var(par, v0)*prePy_0(par, y, y0, v0)*(pdf_c_0(par, y, y0, v2, v0) - 
			pdf_c_0(par, y, y0, v1, v0));
}

inline double Integ2_0(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	return Integ0_0(par, y, y0, v1, v2, v0)*cond_v_var(par, v0) + Integ1_0(par, y, y0, v1, v2, v0)*cond_v_mean_0(par, y, y0, v0) - 
		prePy_0(par, y, y0, v0)*cond_v_var(par, v0)*(pdf_c_0(par, y, y0, v2, v0)*v2 - pdf_c_0(par, y, y0, v1, v0)*v1);
}

inline double Integ3_0(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	return 2*Integ1_0(par, y, y0, v1, v2, v0)*cond_v_var(par, v0) + 
		Integ2_0(par, y, y0, v1, v2, v0)*cond_v_mean_0(par, y, y0, v0) - 
		prePy_0(par, y, y0, v0)*cond_v_var(par, v0)*(pdf_c_0(par, y, y0, v2, v0)*v2*v2 - 
			pdf_c_0(par, y, y0, v1, v0)*v1*v1);
}

inline double Integ0_n(const stdVec& par, double y, double y0, double v1, double v2, double v0, int n) {
	return prePy_n(par, y, y0, v0, n) * (cdf_c_n(par, y, y0, v2, v0, n) - cdf_c_n(par, y, y0, v1, v0, n));
}

inline double Integ1_n(const stdVec& par, double y, double y0, double v1, double v2, double v0, int n) {
	return Integ0_n(par, y, y0, v1, v2, v0, n)*cond_v_mean_n(par, y, y0, v0, n)  - 
		cond_v_var(par, v0)*prePy_n(par, y, y0, v0, n)*(pdf_c_n(par, y, y0, v2, v0, n) - 
			pdf_c_n(par, y, y0, v1, v0, n));
}

inline double Integ2_n(const stdVec& par, double y, double y0, double v1, double v2, double v0, int n) {
	return Integ0_n(par, y, y0, v1, v2, v0, n)*cond_v_var(par, v0) + Integ1_n(par, y, y0, v1, v2, v0, n)*cond_v_mean_n(par, y, y0, v0, n) - 
		prePy_n(par, y, y0, v0, n)*cond_v_var(par, v0)*(pdf_c_n(par, y, y0, v2, v0, n)*v2 - pdf_c_n(par, y, y0, v1, v0, n)*v1);
}

inline double Integ3_n(const stdVec& par, double y, double y0, double v1, double v2, double v0, int n) {
	return 2*Integ1_n(par, y, y0, v1, v2, v0, n)*cond_v_var(par, v0) + 
		Integ2_n(par, y, y0, v1, v2, v0, n)*cond_v_mean_n(par, y, y0, v0, n) - 
		prePy_n(par, y, y0, v0, n)*cond_v_var(par, v0)*(pdf_c_n(par, y, y0, v2, v0, n)*v2*v2 - 
			pdf_c_n(par, y, y0, v1, v0, n)*v1*v1);
}

inline double Integ0(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double ret = Integ0_0(par, y, y0, v1, v2, v0);
	double fact = par[5]*Delta;
	for(int i = 1; i <= njumps; ++i) {
		ret += Integ0_n(par, y, y0, v1, v2, v0, i) * fact;
		fact *= par[5]*Delta/i;
	}
	return exp(-par[5]*Delta) * ret;
}

inline double Integ1(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double ret = Integ1_0(par, y, y0, v1, v2, v0);
	double fact = par[5]*Delta;
	for(int i = 1; i <= njumps; ++i) {
		ret += Integ1_n(par, y, y0, v1, v2, v0, i) * fact;
		fact *= par[5]*Delta/i;
	}
	return exp(-par[5]*Delta) * ret;
}

inline double Integ2(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double ret = Integ2_0(par, y, y0, v1, v2, v0);
	double fact = par[5]*Delta;
	for(int i = 1; i <= njumps; ++i) {
		ret += Integ2_n(par, y, y0, v1, v2, v0, i) * fact;
		fact *= par[5]*Delta/i;
	}
	return exp(-par[5]*Delta) * ret;
}

inline double Integ3(const stdVec& par, double y, double y0, double v1, double v2, double v0) {
	double ret = Integ3_0(par, y, y0, v1, v2, v0);
	double fact = par[5]*Delta;
	for(int i = 1; i <= njumps; ++i) {
		ret += Integ3_n(par, y, y0, v1, v2, v0, i) * fact;
		fact *= par[5]*Delta/i;
	}
	return exp(-par[5]*Delta) * ret;
}

inline double prePy(const stdVec& par, double y, double y0, double v0) {
	double ret = prePy_0(par, y, y0, v0);
	double fact = par[5]*Delta; 
	for(int i = 1; i <= njumps; ++i) {
		ret += prePy_n(par, y, y0, v0, i) * fact;
		fact *= par[5]*Delta/i;
	}
	return exp(-par[5]*Delta) * ret;
}

// lagrange cubic interpolation, coefficients stored in coef
template< typename F >
void interCoef(F& f, double x, double interval, stdVec& coef) {
	double f1 = f(x);
	double f2 = f(x+ interval);
	double f3 = f(x+ 2*interval);
	double f4 = f(x+ 3*interval);

	coef[0] = f1 + ((11*f1 - 18*f2 + 9*f3 - 
		2*f4)*x)/(6.*interval) + ((f1 - (5*f2)/2. + 2*f3 - 
		f4/2.)*pow(x,2))/pow(interval,2) + ((f1 - 3*f2 + 
		3*f3 - f4)*pow(x,3))/(6.*pow(interval,3));

	coef[1] = (-11*f1 + 18*f2 - 9*f3 +
	 2*f4)/(6.*interval) + ((-2*f1 + 5*f2 - 
	 	4*f3 + f4)*x)/pow(interval,2) + ((-f1 + 3*f2 - 
	 	3*f3 + f4)*pow(x,2))/(2.*pow(interval,3));

	coef[2] = (f1 - (5*f2)/2. + 2*f3 - 
		f4/2.)/pow(interval,2) + ((f1 - 3*f2 + 3*f3 - 
			f4)*x)/(2.*pow(interval,3));
	
	coef[3] = (-f1 + 3*f2 - 3*f3 + f4)/(6.*pow(interval,3));
}


// find the nodes for interpolation
// given the (nbasis+1) points [-0.8. 0.4]
// generate from stationary gamma distribution
void findNodes(const stdVec& par, stdVec& interNodes) {
	// create the (nbasis+1) nodes we want to interpolate
	interNodes = myNodes;
	for(int i = 0; i < nbasis+1; ++i) 
		interNodes[i] = exp(log(par[2]) - interNodes[i]*log(sqrt(par[2]*par[3]*par[3]/(2*par[1])))); 	
}

/*
void t1(stdMat& beta, const stdVec& interNodes, int k) {
	for(int j = max(k-3, 0); j < min(k+3, nbasis); ++j) {
			stdVec tmp(4);
			double stepSize = (interNodes[j+1] - interNodes[j])/3.0;
			
			interCoef(I0, interNodes[j], stepSize, tmp);
			beta[k][j] = tmp[0];          beta[k][j+nbasis] = tmp[1];
			beta[k][j+2*nbasis] = tmp[2]; beta[k][j+3*nbasis] = tmp[3];

			interCoef(I1, interNodes[j], stepSize, tmp);
			beta[k+nbasis][j] = tmp[0];          beta[k+nbasis][j+nbasis] = tmp[1];
			beta[k+nbasis][j+2*nbasis] = tmp[2]; beta[k+nbasis][j+3*nbasis] = tmp[3];
		}
}

void t2(stdMat& beta, const stdVec& interNodes, int k) {
	for(int j = max(k-3, 0); j < min(k+3, nbasis); ++j) {
			stdVec tmp(4);
			double stepSize = (interNodes[j+1] - interNodes[j])/3.0;

			interCoef(I2, interNodes[j], stepSize, tmp);
			beta[k+2*nbasis][j] = tmp[0];          beta[k+2*nbasis][j+nbasis] = tmp[1];
			beta[k+2*nbasis][j+2*nbasis] = tmp[2]; beta[k+2*nbasis][j+3*nbasis] = tmp[3];

			interCoef(I3, interNodes[j], stepSize, tmp);
			beta[k+3*nbasis][j] = tmp[0];          beta[k+3*nbasis][j+nbasis] = tmp[1];
			beta[k+3*nbasis][j+2*nbasis] = tmp[2]; beta[k+3*nbasis][j+3*nbasis] = tmp[3];
		}
}
*/


// find the coef vector alpha, expansion of marginal transition density
// and the coef matrix beta, expansion of B_k
void BasisCoef(stdVec& alpha, stdMat& beta, double y, double y0, const stdVec& par, const stdVec& interNodes) {
	// calculate alpha first
	auto p1_y = bind(prePy, par, y, y0, placeholders::_1); // the marginal density functional object, i.e. prestfun

	//#pragma omp parallel shared(interNodes, alpha, beta, par) 
	//#pragma omp parallel for 
	#pragma omp parallel for default(shared) num_threads(1) 
	for(int k = 0; k < nbasis; ++k) {
		stdVec tmp(4);
		double stepSize = (interNodes[k+1] - interNodes[k])/3.0;
		interCoef(p1_y, interNodes[k], stepSize, tmp);
		alpha[k] = tmp[0];
		alpha[k+nbasis] = tmp[1];
		alpha[k+2*nbasis] = tmp[2];
		alpha[k+3*nbasis] = tmp[3];

		// initialize the four integral
		// then interpolate them
		auto I0 = bind(Integ0, par, y, y0, interNodes[k], interNodes[k+1], placeholders::_1);
		auto I1 = bind(Integ1, par, y, y0, interNodes[k], interNodes[k+1], placeholders::_1);
		auto I2 = bind(Integ2, par, y, y0, interNodes[k], interNodes[k+1], placeholders::_1);
		auto I3 = bind(Integ3, par, y, y0, interNodes[k], interNodes[k+1], placeholders::_1);
		// abs(k-j) > 3 are dropped because they're very small in abs
		for(int j = max(k-3, 0); j < min(k+3, nbasis); ++j) {
			// equals to 
			//if(abs(k-j) > 3)
			//	continue;
			stepSize = (interNodes[j+1] - interNodes[j])/3.0;
			
			interCoef(I0, interNodes[j], stepSize, tmp);
			beta[k][j] = tmp[0];          beta[k][j+nbasis] = tmp[1];
			beta[k][j+2*nbasis] = tmp[2]; beta[k][j+3*nbasis] = tmp[3];

			interCoef(I1, interNodes[j], stepSize, tmp);
			beta[k+nbasis][j] = tmp[0];          beta[k+nbasis][j+nbasis] = tmp[1];
			beta[k+nbasis][j+2*nbasis] = tmp[2]; beta[k+nbasis][j+3*nbasis] = tmp[3];

			interCoef(I2, interNodes[j], stepSize, tmp);
			beta[k+2*nbasis][j] = tmp[0];          beta[k+2*nbasis][j+nbasis] = tmp[1];
			beta[k+2*nbasis][j+2*nbasis] = tmp[2]; beta[k+2*nbasis][j+3*nbasis] = tmp[3];

			interCoef(I3, interNodes[j], stepSize, tmp);
			beta[k+3*nbasis][j] = tmp[0];          beta[k+3*nbasis][j+nbasis] = tmp[1];
			beta[k+3*nbasis][j+2*nbasis] = tmp[2]; beta[k+3*nbasis][j+3*nbasis] = tmp[3];
		}
	} 
}

// calculate the initial moments

void initialMoments(const stdVec& par, stdVec& initialMoments, const stdVec& nodes) {
	// create vector to hold the scaled interval nodes
	stdVec interNodes(nbasis+1);
	// scale and shape parameter of gamma distribution
	double beta = 0.5*par[3]*par[3]/par[1];
	double shape = 2*par[1]*par[2]/pow(par[3], 2);
	// scale the nodes by 1/scale
	transform(nodes.begin(), nodes.end(), interNodes.begin(), [=](double x) { return x/beta; });
	//#pragma omp parallel shared(initialMoments, shape, interNodes, beta)
	#pragma omp parallel for default(shared) num_threads(1)
	for(int i = 0; i < nbasis; ++i) {
		initialMoments[i] = gamma_p(shape, interNodes[i+1]) - gamma_p(shape, interNodes[i]);
		initialMoments[nbasis+i] = shape*beta*(gamma_p(shape+1, interNodes[i+1]) - gamma_p(shape+1, interNodes[i]));
		initialMoments[2*nbasis+i] = shape*(shape+1)*beta*beta * (gamma_p(shape+2, interNodes[i+1]) - gamma_p(shape+2, interNodes[i]));
		initialMoments[3*nbasis+i] = shape*(shape+1)*(shape+2)*beta*beta*beta * (gamma_p(shape+3, interNodes[i+1]) - gamma_p(shape+3, interNodes[i]));
	}
}

double nll(const stdVec& par, const stdVec& y) {
	int N = y.size() - 1;
	// M_{k, i}^{(l)}
	stdMat filterMoment(N+1, stdVec(4*nbasis)); 
	
	stdVec initm(4*nbasis);
	stdVec nodes(nbasis+1);
	findNodes(par, nodes);
	initialMoments(par, initm, nodes);
	filterMoment[0] = initm;
	
	stdVec Li(N);
	stdVec alpha(4*nbasis);
	stdMat beta(4*nbasis, stdVec(4*nbasis));

	int i, j;
	double ret;
	double count_error = 0;
	for(i = 0; i < N; ++i) {
		BasisCoef(alpha, beta, y[i+1], y[i], par, nodes);
		Li[i] = inner_product(filterMoment[i].begin(), filterMoment[i].end(), alpha.begin(), 0.0);
		ret -= log(Li[i]);
		
		for(j = 0; j < 4*nbasis; ++j) { 
			filterMoment[i+1][j] = inner_product(beta[j].begin(), beta[j].end(), filterMoment[i].begin(), 0.0) / Li[i];
		}
	}
	return ret;
}

int main() {
	stdVec par = {0.05, 2, 0.2, 0.25, -0.8, 0, 0, 0};
	stdVec y = {4.6051702,4.6194641,4.5787097,4.5914551,4.6029164,4.5893197,4.5819875,4.606393,4.6063427,4.5909786,4.6030319,4.6037265,4.6085072,4.6083776,4.6098088,4.6184302,4.5884666,4.5866887,4.5745194,4.5804907,4.5991651,4.6296578,4.6185303,4.6292363,4.6223832,4.6140163,4.5932847,4.5940766,4.6392253,4.6402955,4.6246525,4.6329189,4.6567086,4.6436899,4.6463981,4.6545833,4.6419823,4.6427593,4.6205012,4.6621959,4.6614725,4.6971056,4.6813843,4.6613061,4.6280359,4.6255658,4.5843473,4.602707,4.6289575,4.6179757,4.6308881,4.6188019,4.6192763,4.6032201,4.6017877,4.6181296,4.6087955,4.594717,4.5964834,4.6115526,4.5750188,4.5864394,4.5959735,4.6057777,4.6073259,4.6020913,4.5985116,4.6157218,4.5865242,4.5850392,4.5714552,4.5653768,4.5689057,4.579034,4.5926883,4.5797717,4.5738295,4.5561819,4.5999294,4.612857,4.6286063,4.6310328,4.6351877,4.6327603,4.6572814,4.6415582,4.611428,4.61694,4.6244043,4.6428807,4.6372718,4.6483665,4.676614,4.682402,4.6743424,4.6623364,4.6724773,4.6675933,4.6729259,4.7032723,4.6877363,4.7069819,4.7061388,4.6904546,4.6756661,4.6633886,4.6915157,4.6805667,4.6799109,4.6924511,4.6700996,4.6605403,4.6784427,4.6845075,4.6827808,4.687037,4.653718,4.6304225,4.6326216,4.6548115,4.6361641,4.6060547,4.5576224,4.5512473,4.5873674,4.5790121,4.5850632,4.5941382,4.5829025,4.621832,4.5748235,4.6018109,4.578725,4.6436967,4.6678826,4.6651678,4.65789,4.6825284,4.7100294,4.688825,4.6865675,4.697394,4.6898871,4.6903812,4.6726531,4.6803761,4.7244369,4.7448464,4.7542344,4.7305176,4.7189503,4.7526131,4.7356572,4.7559062,4.7677838,4.7793433,4.8134381,4.8066634,4.8069393,4.7958663,4.8066511,4.8083684,4.803281,4.8021195,4.8187537,4.8122124,4.7911116,4.7995065,4.7996387,4.8114697,4.8210114,4.8325883,4.8343761,4.8380987,4.8475515,4.8520925,4.8578463,4.8688549,4.8560987,4.8580085,4.8869389,4.8646037,4.8563215,4.8511439,4.8479998,4.8450629,4.7824271,4.7787975,4.7881921,4.7941404,4.7752344,4.771889,4.763198,4.7680585,4.7586063,4.7641828,4.7618846,4.7416919,4.7364424,4.7442576,4.7216649,4.7324746,4.7333858,4.7443262,4.7466382,4.7273278,4.7296776,4.7169459,4.7146317,4.6955679,4.6871827,4.6643963,4.6560757,4.6162386,4.6528751,4.6778241,4.6655508,4.6523493,4.6220536,4.5825029,4.5815598,4.5715193,4.5725374,4.5689871,4.5689509,4.5509075,4.5683878,4.5769302,4.5636044,4.5926461,4.5783941,4.5924038,4.5904616,4.5923116,4.5950425,4.5886458,4.6089715,4.6286693,4.6396439,4.6180222,4.616155,4.6453947,4.6604664,4.654596,4.6241795,4.6372185,4.5971862,4.6102855,4.6149502,4.6331067,4.629256,4.6427416,4.630131,4.6334021,4.6068088,4.6178391,4.5919613,4.6093773,4.613289,4.5996805,4.5712488,4.5694448,4.5852463,4.5599186,4.5884485,4.5759153,4.5903516,4.5801071,4.5766107,4.580001,4.5637865,4.5749192,4.5568814,4.5987294,4.6015087,4.5896502,4.5879269,4.596172,4.5656115,4.5217264,4.5273637,4.4639603,4.4698486,4.461827,4.4611853,4.4584465,4.453696,4.481514,4.4822096,4.4893567,4.4691896,4.4861856,4.5311377,4.528875,4.5470011,4.576984,4.5713188,4.583144,4.6083433,4.5962519,4.6083188,4.5979763,4.595435,4.6258694,4.6033688,4.5733307,4.5339653,4.5588774,4.6056068,4.6040921,4.5903765,4.5718122,4.5380431,4.5629959,4.5473108,4.5256659,4.5385891,4.5554792,4.5908568,4.609427,4.6011671,4.5926818,4.6148571,4.6381891,4.6327459,4.6237516,4.5549864,4.5432124,4.5169961,4.5342231,4.4923996,4.5109156,4.5180595,4.5316216,4.5378072,4.5873986,4.5962878,4.6472854,4.66028,4.6329751,4.6252394,4.6111615,4.5840506,4.5831034,4.5834639,4.5838566,4.5928528,4.5881949,4.5792229,4.5682386,4.5664462,4.603161,4.6206667,4.665193,4.6466976,4.6373279,4.665023,4.7005366,4.6906451,4.699755,4.7145446,4.7437969,4.7753036,4.7573647,4.7438775,4.7379958,4.7249774,4.7427865,4.7052059,4.734246,4.7534865,4.7540193,4.7456576,4.7439513,4.7036951,4.7059426,4.7031775,4.7174513,4.7536968,4.7603392,4.7912955,4.7946151,4.7959046,4.8241147,4.8096515,4.8230817,4.8453084,4.8648939,4.8780866,4.8756031,4.8826237,4.8890248,4.8739512,4.8618208,4.8768171,4.8894702,4.9043219,4.9138032,4.9147879,4.9059049,4.9101275,4.8826479,4.8979647,4.9216745,4.9383237,4.924966,4.9288373,4.9462817,4.9335897,4.9776768,4.9858238,4.9795429,4.9876294,4.9851871,4.9869858,4.9716336,4.9756202,4.9277457,4.8947933,4.9214173,4.9126867,4.9192642,4.9275525,4.9245778,4.9377647,4.9315301,4.9259692,4.9419096,4.9440819,4.9511659,4.9530356,4.946606,4.9631917,4.9776839,4.9422193,4.9667153,4.984542,4.9637515,4.9667102,4.9771416,4.9809825,4.9864532,4.9907897,4.9856967,4.9614919,4.9575643,4.9364645,4.9225821,4.900139,4.911873,4.8878871,4.8929158,4.9115893,4.9127346,4.928334,4.9190766,4.9230175,4.9186999,4.9040949,4.8872712,4.882465,4.9179407,4.8967923,4.8655938,4.8793122,4.86428,4.9059463,4.8742644,4.8903043,4.8837126,4.8649619,4.8559992,4.836795,4.8141799,4.7996406,4.8112392,4.8056191,4.8042521,4.8174445,4.8015287,4.7892547,4.7977896,4.796562,4.7880285,4.7842701,4.7960562,4.7912628,4.8110237,4.8156595,4.8108652,4.79865,4.8028632,4.7963571,4.7974164,4.7822524,4.7503857,4.7606413,4.7514481,4.7420287,4.7646363,4.787242,4.798284,4.7821849,4.7755627,4.7745378,4.7649964,4.795397,4.8004385,4.782713,4.7217447,4.7302103,4.6963617,4.710126,4.695103,4.7207001,4.7079461,4.7522597,4.7407479,4.7639724,4.7370624,4.7372148,4.737284,4.7301024,4.7146394,4.7039707,4.6946433,4.6986367,4.7089802,4.6845374,4.6978529,4.7129154,4.6833241,4.6687174,4.659759,4.606498,4.5957375,4.6152682,4.6522671,4.6336721,4.6034175,4.5968781,4.5874318,4.5676588,4.5451705,4.52754,4.5401965,4.5157659,4.5428107,4.5638203,4.5819751,4.6062525,4.6138168,4.6207399,4.599723,4.601061,4.579872,4.5716347,4.5762724,4.6399431,4.6547033,4.6410651,4.625379,4.6147978,4.5792103,4.546107,4.5759663,4.5368003,4.5129426,4.5239124,4.5179202,4.4969486,4.5124079,4.5325054,4.5702134,4.5341434,4.5193732,4.4989244,4.5310655,4.5328254,4.5574179,4.5493274,4.5500742,4.5666593,4.5651267,4.5823548,4.5700071,4.5457455,4.5416496,4.5574268,4.5689435,4.5762408,4.5676123,4.5660672,4.5729076,4.554578,4.5418883,4.5540897,4.5634966,4.540515,4.5825567,4.5834567,4.6112918,4.590335,4.6010727,4.6047109,4.6186807,4.618525,4.6193143,4.6289614,4.6652422,4.6710094,4.7009136,4.7095696,4.6916312,4.7139176,4.7414035,4.7460021,4.752234,4.7810921,4.770777,4.7703062,4.7539316,4.7410799,4.752793,4.7687927,4.7315724,4.7120843,4.6964266,4.7094517,4.6898499,4.6791256,4.6745462,4.6705396,4.6667178,4.6720862,4.6786179,4.6886303,4.6863853,4.6897097,4.7139632,4.7518528,4.7831462,4.7574865,4.7635261,4.7309415,4.7541399,4.7594343,4.7790945,4.7594892,4.7548251,4.7594376,4.7651037,4.7299061,4.7683355,4.744659,4.6949978,4.7208114,4.6983179,4.7071216,4.6955805,4.707077,4.7265474,4.7175933,4.7090903,4.6936715,4.6891581,4.6801617,4.6677705,4.6558048,4.6452548,4.6359709,4.6229895,4.6014985,4.6304369,4.6118363,4.647948,4.6394022,4.6307368,4.6257131,4.6391599,4.6183246,4.6259192,4.6168042,4.640857,4.6593427,4.6179985,4.6154036,4.5618966,4.5584336,4.4887855,4.4723225,4.4753267,4.4617298,4.5047334,4.5054586,4.5094351,4.5360467,4.5298929,4.5257925,4.5096632,4.4977791,4.489972,4.4942688,4.4962713,4.5078164,4.494082,4.4693401,4.4589778,4.4485301,4.4423303,4.4241136,4.4559264,4.4565181,4.4759706,4.4766985,4.5097901,4.4954974,4.5310616,4.5089746,4.5379477,4.5180824,4.5230913,4.5539828,4.5585471,4.5739851,4.5734957,4.5534159,4.5473469,4.527444,4.5228993,4.5360902,4.5405669,4.5764136,4.5583187,4.5699083,4.5640819,4.5667332,4.5940227,4.6070909,4.6247508,4.6390189,4.6647542,4.645922,4.6457896,4.6461203,4.6267934,4.6308189,4.6607344,4.6609122,4.661729,4.6697597,4.6520194,4.6461815,4.6272573,4.6038127,4.6123438,4.6141644,4.595511,4.6163092,4.6571925,4.6231646,4.6269612,4.6342422,4.631403,4.6355698,4.6034913,4.5910268,4.6158111,4.6171665,4.5998509,4.5724408,4.6006387,4.5913595,4.538745,4.5547068,4.5246997,4.5172537,4.5024283,4.4980871,4.4738406,4.4998888,4.4846938,4.4482874,4.4603833,4.4556242,4.4761828,4.4777264,4.5146264,4.5306975,4.5613731,4.559621,4.5545378,4.5447109,4.5464573,4.5441603,4.5436663,4.5634321,4.5768053,4.5867152,4.5623202,4.5408332,4.5724499,4.5845673,4.5648885,4.5762892,4.5474812,4.5435066,4.561056,4.5549741,4.5247373,4.5015798,4.5059702,4.5093579,4.5095982,4.487031,4.5138253,4.4823969,4.4474713,4.4643606,4.5182834,4.5248364,4.5460335,4.5685355,4.5679031,4.5402718,4.5523302,4.538588,4.5225384,4.5179364,4.5252961,4.5392584,4.5322755,4.510796,4.4870446,4.4972151,4.5223543,4.5202217,4.5412951,4.5170425,4.5245741,4.5249188,4.5251977,4.5173794,4.5129006,4.4799532,4.4916571,4.4644138,4.4693419,4.4443623,4.4317069,4.4374879,4.4208931,4.423371,4.4192748,4.3864674,4.4080651,4.4152905,4.4678805,4.4551206,4.4624119,4.4571176,4.485231,4.4722727,4.4684085,4.5006987,4.4947801,4.4950878,4.4646967,4.5027085,4.5297118,4.519615,4.5380361,4.5154796,4.5317843,4.544093,4.5543513,4.5344455,4.5447643,4.5346609,4.5190012,4.524988,4.531669,4.526062,4.5446921,4.5101369,4.5167331,4.5057947,4.4962894,4.4994983,4.4982584,4.5102891,4.5200043,4.5259684,4.5209106,4.5169985,4.4940264,4.5036151,4.4785031,4.4944344,4.4831479,4.5036082,4.5215677,4.4672697,4.4528827,4.4556639,4.4500127,4.3996919,4.3950984,4.4294825,4.4135852,4.4127598,4.4286344,4.44829,4.4556077,4.4525403,4.4284404,4.4361459,4.4179887,4.4110205,4.3877752,4.4096273,4.4256374,4.4389132,4.4291875,4.4599681,4.4282805,4.4233363,4.4368194,4.4432215,4.4616522,4.4659195,4.4488284,4.4270339,4.4536924,4.4828236,4.4819571,4.447816,4.4143021,4.414082,4.4006659,4.369653,4.4301341,4.4527965,4.466124,4.4151707,4.456952,4.4693467,4.432147,4.46223,4.4296901,4.4200832,4.4416553,4.470955,4.4649036,4.4686929,4.4927279,4.4790549,4.499243,4.5058706,4.5353859,4.5624506,4.5609305,4.5408159,4.5594009,4.5576369,4.5784528,4.5749793,4.5770438,4.588135,4.615418,4.6003581,4.6022932,4.5969117,4.6210469,4.5954908,4.6092676,4.5740613,4.5660551,4.5773584,4.5795843,4.5805867,4.6061281,4.6108213,4.6221166,4.6142314,4.6117394,4.6378066,4.6473028,4.6517075,4.646456,4.6647881,4.6500794,4.6391091,4.6455952,4.6575058,4.646459,4.6336136,4.6296883,4.6498705,4.6305003,4.6399666};
	double res = 0;
	int times = 1;
	clock_t start = clock();
	for(size_t i = 0; i < times; ++i) {
		res = nll(par, y);
	}
	clock_t end = clock();
	double seconds = (end - start) / (double) CLOCKS_PER_SEC;

	cout << setprecision(12) << setw(12);
	cout << res << endl;
	cout << "time cost: " << seconds / (double) times << endl;
}