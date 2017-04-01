%% value of cash flow pe ratio
%cum dividend price

% the cash flow growth: mu_e + phi_e*(x_t - x_bar)+sig*sig_e *eta_c-
% 0.5*sig^2*sig_e^2*sig_c^2 +sig*eta_id - 0.5*sig^2*sig_eid^2+(d_e*J - (exp(d_e)-1)*lambda*x_t)

pe = zeros(nx, 1);

pe0 = pe;

pe1 = pe0;

diff = 1;
while diff > sqrt(eps)    
%     vexp = v0*px';
%     beta = inv(X'*X)*X'*v0;

    consG = exp(mu_e + phi_e*(x - x_bar)-0.5*sig.^2*sig_e^2*sig_c^2 ...
    +0.5*(sig_e-gamma)^2*sig.^2*sig_c^2);
 
    % px.*M risk neutral prob
    pe1 = 1-tau+consG.*((px.*M)*pe0);         

    diff = norm(pe1(:) - pe0(:));
    
    %disp(diff);
    
    pe0 = pe1;
    
end

pe = pe0;