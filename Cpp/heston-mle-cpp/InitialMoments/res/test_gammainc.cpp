#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <boost/math/special_functions/gamma.hpp>

using boost::math;

#include "Parameter.hpp"

using namespace std;

typedef vector< double > stdVec;

int main() {
    double mu = 0.05;
    double kappa = 2;
    double theta = 0.2;
    double xi = 0.25;
    double rho = -0.8;

    cout << "test incomplete gamma lower function in boost: " << endl;
	double z1 = 2*kappa*theta/(xi*xi);
	double a1 = 0.154221*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.375);
	double z2 = 2*kappa*theta/(xi*xi);
	double a2 = 0.151572*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.4);
    cout << std::setprecision(16) << boost::math::gamma_p(z1, a1) - boost::math::gamma_p(z2, a2) << endl;
	return 0;
}


/*
	stdVec initialMoments(82); 4 x 21
*/ 

stdVec initialMoments(Paramater Par) {
	mu = Par.mu;
	kappa = Par.kappa;
	theta = Par.theta;
	xi = Par.xi;
	rho = Par.rho;

	initialMoments[0] = gamma_p(2*kappa*theta/(xi*xi), 0.154221*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.375)) - 
		gamma_p(2*kappa*theta/(xi*xi), 0.151572*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.4));
	initialMoments[1] = gamma_p(2*kappa*theta/(xi*xi), 0.156917*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.35)) - 
		gamma_p(2*kappa*theta/(xi*xi), 0.154221*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.375));
	initialMoments[2] = gamma_p(2*kappa*theta/(xi*xi), 0.15966*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.325)) - 
		gamma_p(2*kappa*theta/(xi*xi), 0.151572*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.4));
	initialMoments[3] = gamma_p(2*kappa*theta/(xi*xi), 0.156917*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.3)) - 
		gamma_p(2*kappa*theta/(xi*xi), 0.154221*10*kappa*theta/(xi*xi)*pow(1/kappa*theta*xi*xi, 0.375));
		

	0.154221.^0.375E0        0.151572.^0.4E0);
	0.156917.^0.35E0        0.154221.^0.375E0);
	0.15966.^0.325E0        0.156917.^0.35E0);
	0.16245.^0.3E0        0.15966.^0.325E0);
	0.16529.^0.275E0        0.16245.^0.3E0);
	0.168179.^0.25E0        0.16529.^0.275E0);
	0.171119.^0.225E0        0.168179.^0.25E0);
	0.17411.^0.2E0        0.171119.^0.225E0);
	0.177154.^0.175E0        0.17411.^0.2E0);
	0.18025.^0.15E0        0.177154.^0.175E0);
	0.183401.^0.125E0        0.18025.^0.15E0);
	0.186607.^0.1E0        0.183401.^0.125E0);
	0.189868.^0.75E-1        0.186607.^0.1E0);
	0.193187.^0.5E-1        0.189868.^0.75E-1);
	0.196564.^0.25E-1        0.193187.^0.5E-1);
	2*kappa*theta/(xi*xi)        0.196564.^0.25E-1);
	initialMoments[17]=(-1)*GammahRegularized(2*kappa*theta/(xi*xi),0,2*kappa*theta/(xi*xi))+GammahRegularized(2*kappa*theta/(xi*xi),0,0.203496.^(-0.25E-1));
	0.207053.^(-0.5E-1)        0.203496.^(-0.25E-1));
	0.214355.^(-0.1E0)        0.207053.^(-0.5E-1));
	0.221914.^(-0.15E0)        0.214355.^(-0.1E0));
	0.22974.^(-0.2E0)        0.221914.^(-0.15E0));
	



















	initialMoments[22]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.154221.^0.375E0)+Gammah(1+2*kappa*theta/(xi*xi),0.151572.^0.4E0));
	initialMoments[23]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.156917.^0.35E0)+Gammah(1+2*kappa*theta/(xi*xi),0.154221.^0.375E0));
	initialMoments[24]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.15966.^0.325E0)+Gammah(1+2*kappa*theta/(xi*xi),0.156917.^0.35E0));
	initialMoments[25]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.16245.^0.3E0)+Gammah(1+2*kappa*theta/(xi*xi),0.15966.^0.325E0));
	initialMoments[26]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.16529.^0.275E0)+Gammah(1+2*kappa*theta/(xi*xi),0.16245.^0.3E0));
	initialMoments[27]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.168179.^0.25E0)+Gammah(1+2*kappa*theta/(xi*xi),0.16529.^0.275E0));
	initialMoments[28]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.171119.^0.225E0)+Gammah(1+2*kappa*theta/(xi*xi),0.168179.^0.25E0));
	initialMoments[29]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.17411.^0.2E0)+Gammah(1+2*kappa*theta/(xi*xi),0.171119.^0.225E0));
	initialMoments[30]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.177154.^0.175E0)+Gammah(1+2*kappa*theta/(xi*xi),0.17411.^0.2E0));
	initialMoments[31]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.18025.^0.15E0)+Gammah(1+2*kappa*theta/(xi*xi),0.177154.^0.175E0));
	initialMoments[32]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.183401.^0.125E0)+Gammah(1+2*kappa*theta/(xi*xi),0.18025.^0.15E0));
	initialMoments[33]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.186607.^0.1E0)+Gammah(1+2*kappa*theta/(xi*xi),0.183401.^0.125E0));
	initialMoments[34]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.189868.^0.75E-1)+Gammah(1+2*kappa*theta/(xi*xi),0.186607.^0.1E0));
	initialMoments[35]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.193187.^0.5E-1)+Gammah(1+2*kappa*theta/(xi*xi),0.189868.^0.75E-1));
	initialMoments[36]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.196564.^0.25E-1)+Gammah(1+2*kappa*theta/(xi*xi),0.193187.^0.5E-1));
	initialMoments[37]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),2*kappa*theta/(xi*xi))+Gammah(1+2*kappa*theta/(xi*xi),0.196564.^0.25E-1));
	initialMoments[38]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*(Gammah(1+2*kappa*theta/(xi*xi),2*kappa*theta/(xi*xi))+(-1)*Gammah(1+2*kappa*theta/(xi*xi),0.203496.^(-0.25E-1)));
	initialMoments[39]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.207053.^(-0.5E-1))+Gammah(1+2*kappa*theta/(xi*xi),0.203496.^(-0.25E-1)));
	initialMoments[40]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.214355.^(-0.1E0))+Gammah(1+2*kappa*theta/(xi*xi),0.207053.^(-0.5E-1)));
	initialMoments[41]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.221914.^(-0.15E0))+Gammah(1+2*kappa*theta/(xi*xi),0.214355.^(-0.1E0)));
	initialMoments[42]=kappa.^(2*kappa*theta/(xi*xi))*theta*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^((-4)*kappa*theta/(xi*xi))*Gammah(1+2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(1+2*kappa*theta/(xi*xi),0.22974.^(-0.2E0))+Gammah(1+2*kappa*theta/(xi*xi),0.221914.^(-0.15E0)));
	initialMoments[43]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.154221.^0.375E0)+Gammah(2+2*kappa*theta/(xi*xi),0.151572.^0.4E0));
	initialMoments[44]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.156917.^0.35E0)+Gammah(2+2*kappa*theta/(xi*xi),0.154221.^0.375E0));
	initialMoments[45]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.15966.^0.325E0)+Gammah(2+2*kappa*theta/(xi*xi),0.156917.^0.35E0));
	initialMoments[46]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.16245.^0.3E0)+Gammah(2+2*kappa*theta/(xi*xi),0.15966.^0.325E0));
	initialMoments[47]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.16529.^0.275E0)+Gammah(2+2*kappa*theta/(xi*xi),0.16245.^0.3E0));
	initialMoments[48]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.168179.^0.25E0)+Gammah(2+2*kappa*theta/(xi*xi),0.16529.^0.275E0));
	initialMoments[49]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.171119.^0.225E0)+Gammah(2+2*kappa*theta/(xi*xi),0.168179.^0.25E0));
	initialMoments[50]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.17411.^0.2E0)+Gammah(2+2*kappa*theta/(xi*xi),0.171119.^0.225E0));
	initialMoments[51]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.177154.^0.175E0)+Gammah(2+2*kappa*theta/(xi*xi),0.17411.^0.2E0));
	initialMoments[52]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.18025.^0.15E0)+Gammah(2+2*kappa*theta/(xi*xi),0.177154.^0.175E0));
	initialMoments[53]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.183401.^0.125E0)+Gammah(2+2*kappa*theta/(xi*xi),0.18025.^0.15E0));
	initialMoments[54]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.186607.^0.1E0)+Gammah(2+2*kappa*theta/(xi*xi),0.183401.^0.125E0));
	initialMoments[55]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.189868.^0.75E-1)+Gammah(2+2*kappa*theta/(xi*xi),0.186607.^0.1E0));
	initialMoments[56]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.193187.^0.5E-1)+Gammah(2+2*kappa*theta/(xi*xi),0.189868.^0.75E-1));
	initialMoments[57]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.196564.^0.25E-1)+Gammah(2+2*kappa*theta/(xi*xi),0.193187.^0.5E-1));
	initialMoments[58]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),2*kappa*theta/(xi*xi))+Gammah(2+2*kappa*theta/(xi*xi),0.196564.^0.25E-1));
	initialMoments[59]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*(Gammah(2+2*kappa*theta/(xi*xi),2*kappa*theta/(xi*xi))+(-1)*Gammah(2+2*kappa*theta/(xi*xi),0.203496.^(-0.25E-1)));
	initialMoments[60]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.207053.^(-0.5E-1))+Gammah(2+2*kappa*theta/(xi*xi),0.203496.^(-0.25E-1)));
	initialMoments[61]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.214355.^(-0.1E0))+Gammah(2+2*kappa*theta/(xi*xi),0.207053.^(-0.5E-1)));
	initialMoments[62]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.221914.^(-0.15E0))+Gammah(2+2*kappa*theta/(xi*xi),0.214355.^(-0.1E0)));
	initialMoments[63]=(1/4)*kappa.^((-2)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(4+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(2+2*kappa*theta/(xi*xi),0.22974.^(-0.2E0))+Gammah(2+2*kappa*theta/(xi*xi),0.221914.^(-0.15E0)));
	initialMoments[64]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.154221.^0.375E0)+Gammah(3+2*kappa*theta/(xi*xi),0.151572.^0.4E0));
	initialMoments[65]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.156917.^0.35E0)+Gammah(3+2*kappa*theta/(xi*xi),0.154221.^0.375E0));
	initialMoments[66]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.15966.^0.325E0)+Gammah(3+2*kappa*theta/(xi*xi),0.156917.^0.35E0));
	initialMoments[67]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.16245.^0.3E0)+Gammah(3+2*kappa*theta/(xi*xi),0.15966.^0.325E0));
	initialMoments[68]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.16529.^0.275E0)+Gammah(3+2*kappa*theta/(xi*xi),0.16245.^0.3E0));
	initialMoments[69]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.168179.^0.25E0)+Gammah(3+2*kappa*theta/(xi*xi),0.16529.^0.275E0));
	initialMoments[70]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.171119.^0.225E0)+Gammah(3+2*kappa*theta/(xi*xi),0.168179.^0.25E0));
	initialMoments[71]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.17411.^0.2E0)+Gammah(3+2*kappa*theta/(xi*xi),0.171119.^0.225E0));
	initialMoments[72]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.177154.^0.175E0)+Gammah(3+2*kappa*theta/(xi*xi),0.17411.^0.2E0));
	initialMoments[73]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.18025.^0.15E0)+Gammah(3+2*kappa*theta/(xi*xi),0.177154.^0.175E0));
	initialMoments[74]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.183401.^0.125E0)+Gammah(3+2*kappa*theta/(xi*xi),0.18025.^0.15E0));
	initialMoments[75]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.186607.^0.1E0)+Gammah(3+2*kappa*theta/(xi*xi),0.183401.^0.125E0));
	initialMoments[76]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.189868.^0.75E-1)+Gammah(3+2*kappa*theta/(xi*xi),0.186607.^0.1E0));
	initialMoments[77]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.193187.^0.5E-1)+Gammah(3+2*kappa*theta/(xi*xi),0.189868.^0.75E-1));
	initialMoments[78]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.196564.^0.25E-1)+Gammah(3+2*kappa*theta/(xi*xi),0.193187.^0.5E-1));
	initialMoments[79]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),2*kappa*theta/(xi*xi))+Gammah(3+2*kappa*theta/(xi*xi),0.196564.^0.25E-1));
	initialMoments[80]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*(Gammah(3+2*kappa*theta/(xi*xi),2*kappa*theta/(xi*xi))+(-1)*Gammah(3+2*kappa*theta/(xi*xi),0.203496.^(-0.25E-1)));
	initialMoments[81]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.207053.^(-0.5E-1))+Gammah(3+2*kappa*theta/(xi*xi),0.203496.^(-0.25E-1)));
	initialMoments[82]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.214355.^(-0.1E0))+Gammah(3+2*kappa*theta/(xi*xi),0.207053.^(-0.5E-1)));
	initialMoments[83]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.221914.^(-0.15E0))+Gammah(3+2*kappa*theta/(xi*xi),0.214355.^(-0.1E0)));
	initialMoments[84]=(1/8)*kappa.^((-3)+2*kappa*theta/(xi*xi))*(kappa/(xi*xi)).^((-2)*kappa*theta/(xi*xi))*xi.^(6+(-4)*kappa*theta/(xi*xi))*Gammah(2*kappa*theta/(xi*xi)).^(-1)*((-1)*Gammah(3+2*kappa*theta/(xi*xi),0.22974.^(-0.2E0))+Gammah(3+2*kappa*theta/(xi*xi),0.221914.^(-0.15E0)));
	
}
