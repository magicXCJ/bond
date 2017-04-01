% monthly parameters calibration
%%% preference parameters:

delta = 0.9982; % discount factor
ksi = 1.5; % IES
gamma = 10; % risk aversion
mu = 0.0015;
sigma_a = 0.0052;
lambda = 0.2;
d = 0.15;
rho = 0.989;
m_w = 0.003;
std_w = 0.001;
nu = m_w^2/std_w^2;
xi = m_w*(1-rho)/nu;
phi = 3;
sigma_d = 0.0052;

par = [delta; ksi; gamma; mu; sigma_a; lambda; d; nu; xi; rho; phi; sigma_d];

[z, zm, A_0c, A_1c, A_0m, A_1m, m0, mw, lambda_w, lambda_e] = modelsolve_poisson(par);
[B0, B1] = bondprice(nu, xi, rho, m0, mw, lambda_w, lambda_e, 60);

k1m = exp(zm)./(1+exp(zm));
k0m = log(1+exp(zm))-k1m.*zm;
rm0 = k0m+(k1m-1)*A_0m - nu*log(1-A_1m*xi)+mu+0.5*phi^2*sigma_a^2+0.5*phi^2*sigma_d^2 +phi*lambda*m_w;
rm1 = (k1m*A_1m-phi*lambda)*rho/(1-(k1m*A_1m-phi*lambda)*xi)-A_1m;
rf0 = -B0(1);
rf1 = -B1(1);
rm0 
rf0
rm1
rf1
rf = (rf0+rf1*m_w)*12
pre = (rm0-rf0+(rm1-rf1)*m_w)*12


