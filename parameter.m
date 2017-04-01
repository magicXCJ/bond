% Quarterly calibration of 1 state model
% simple case no time varying vol
%% preference parameters:
%delta = 0.9745; % discount factor
delta = 0.8;
%delta = 0.997; 
ksi = 1.5; % IES
gamma = 10; % risk aversion
theta = (1-gamma)/(1-1/ksi);

%% consumption growth dynamics:
mu_c = 0.02/4; % unconditional mean of consumption growth
phi_c = -0.2; % loading of consumption on disaster probablity
sig_c = 0.0135; % unconditional standard deviation of shock to c;
d = -0.18;% loss given disaster;
x_bar = 0.0065; %unconditonal mean of disaster probablity
rho = 0.8847; %0.979; % persistence paramter of x 
sig_x = 0.0012; % unconditional standard devidation of shock to x
%lambda =1; % correlation of sigma to x

%firm dynamics
%% solving the firm problem
mu_e = 0.02/4; %drift of earnings growth
%mu_e = 0;
phi_e = -13.6; %loading of earnigs on x
%phi_e = -0.8;
sig_e = 2.5; % systematic vol of earnings growth
%sig_e = 1; 
sig_eid = 5*sig_c; %idiosyncratic vol of earnings growth
d_e = -0.2; %size of disaster of earnings growth
lambda_e = 2; % the probablity of disaster for earings growth take the form lambda*x_t

%% levered equity and bond value;
coupon = 0.02; %coupon rate;
%coupon = 0.4;
tau = 0.14; % coporate tax rate;
%tau = 0;
phi_0 = 0.006; %fixed debt adjustment cost;
%phi_0 =0;
phi_1 = 0.02; %proportional adjustment cost;
%phi_1 =0;
alpha = 0.2; %bankrupcy cost;
%alpha = 0;

kappa = 1;

%% value function parameters
% grids
% aggregate state
nx = 5;
dx = sig_x/sqrt(1-rho^2);
[x,px] = tauchen(nx,x_bar,rho,sig_x,2.5);

% firm choice variable
nb = 20;
%nbp = 10;

maxb = 100;
curvature = 0.5;
% b = linspace(eps, maxb, nb)';
% bp = linspace(eps, maxb, nbp)';
b = linspace(eps^curvature, maxb^curvature, nb)'.^(1/curvature);
%bp = linspace(eps^curvature, maxb^curvature, nbp)'.^(1/curvature);

% To approximate log growth of cash flow
% we need two normal shocks (one systematic, one idiosyncratic and also a poission shock)
% The  poission shock is approximated by p0 = 1-x, p(1) = x;
% turn off idiosyncratic vol 
Nj = 2;
Ne= 5; 
[eta, w] = GaussHermite(Ne); 
eta = eta*sqrt(2);
w = w/sqrt(pi);
Ne_id = 5;
[eta_id, w_id] = GaussHermite(Ne); 
eta_id = eta_id*sqrt(2);
w_id = w_id/sqrt(pi);