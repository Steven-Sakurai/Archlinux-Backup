class Parameter {
    public:
        double mu = 0.05, kappa = 2, theta = 0.2, xi = 0.25, rho = -0.8;
		
		//Parameter() {}
		Parameter(double mu_ = 0.05, double kappa_ = 2, double theta_ = 0.2, double xi_ = 0.25, double rho_ = -0.8):
			mu(mu_), kappa(kappa_), theta(theta_), xi(xi_), rho(rho_) {}
};
