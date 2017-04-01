%% iteration
%log pc ratio!!

pc = zeros(nx, 1);

pc0 = pc;

pc1 = pc0;

% store xk in nz*N
sig = ones(nx,1);

diff = 1;
while diff > sqrt(eps)    

    cons = exp(theta*log(delta)+(1-gamma)*(mu_c+phi_c*(x-x_bar))+(1-gamma)^2*0.5*sig.^2*sig_c^2 ...
    +(exp((1-gamma)*d)-1)*x-(1-gamma)*(exp(d)-1)*x);
    pc1 =log(cons.*(px*(exp(log(exp(pc0)+1)*theta))))/theta;
    
    diff = norm(pc1(:) - pc0(:));
    
    %disp(diff);
    
    pc0 = pc1;
end

pc = pc0;

% the SDF
% consM = exp(theta*log(delta)+(-gamma)*(mu_c+phi_c*(x(t)-x_bar))+(-gamma)*sig*sig_c*eta_c....
%           +(exp((-gamma)*d)-1)*x(t)-(-gamma)*(exp(d)-1)*x(t));
% M = consM*(exp(log(exp(pc(xk))+1)/exp(pc(x)))*(theta-1));

%store the SDF in nz*N matrix as a function of x and xk (don't forget the extra +(-gamma)*sig*sig_c*eta_c term)
M = zeros(nx,nx);
for idx = 1:nx 
       consM = exp(theta*log(delta)+(-gamma)*(mu_c+phi_c*(x(idx)-x_bar))+...
           +(exp((-gamma)*d)-1)*x(idx)-(-gamma)*(exp(d)-1)*x(idx));
       pcg = ((exp(pc)+1)/exp(pc(idx))).^(theta-1);
       M(idx,:) = consM*pcg';    
end
