% shock, state t, state t-1, b_t/e_(t-1)
% 
pb0 = pb_kappa1;
pb1 = pb0;

pv0 = pv_kappa1;
pv1 = pv0;

kappa_iter = linspace(0.99, 0.05, 90);
for iter = 1:90
    kappa = kappa_iter(iter);
    % store iteration constants conditional on kappa
    issue = kron(ones(nzxxb,1),b')-(1-kappa)*kron(ones(1,nb),zxxbgrid(:,6)./kron(ones(nb,1),ge0));
    ne_penalty = (issue<0).*issue*1000000;
    %inertia : the smallest positive bond issurance, fixed issurance cost = 0
    [tmp] = min(abs(issue).*(issue>0),[],2);
    [~,inertia] = min(abs(issue - kron(ones(1,nb),tmp)),[],2);
    %inertia policy function
    inertia_pol = reshape(inertia,[nzxx,nb]);
    
    adjcost = phi_0*(issue>kron(ones(1,nb),tmp))+phi_1*issue;
    div_cons = (1-tau) - ((1-tau)*coupon+kappa)*kron(ones(1,nb),zxxbgrid(:,6)./kron(ones(nb,1),ge0))-adjcost;                
    div_cons = div_cons+ne_penalty;
    
    
    diff_b = 1;
    diff_v = 1;
    diff = 1;
    loop=0;
    
while diff > sqrt(eps)%0.0001
    
    loop = loop+1;
    div = div_cons + kron(ones(nx*nb,1),kron(pb0,ones(nz,1))).*issue;
    expv = pzxx_q * (kron(ge0,ones(1,nb)).*reshape(pv0,nzxx,nb));
    [vmax, polb] = max(div+kron(ones(nb,1),expv),[],2);
    pv1 = max(0,vmax); 

    diff_v = norm(pv1(:) - pv0(:));
    %diff_v = sumabs(pv1(:) - pv0(:));
    %disp(diff_v);
    pv0 = pv1;

    polbreal = b(polb);
% given default policy, calculate the bond price, policy consistent?
%disp(pv0)
    pv0_r = reshape(pv0,nzxx,nb);
    
    pb0_expand = kron(ones(nx,1),kron(pb0,ones(nz,1)));
    
    pb_inertia = zeros(nzxx,nb);
    for i = 1:nzxx;
        tmp = pb0_expand(i,:);
        pb_inertia(i,:) = tmp(inertia_pol(i,:));
    end

    Id = (pv0_r<=0);
    ndef = (1-Id).*(coupon+kappa+(1-kappa)*pb_inertia);
    def = Id.*defpayoff;
    pb1 = px_zxx_q*(ndef+def);
    diff_b = norm(pb1(:) - pb0(:));
%     disp('the bond function error is')
%     disp(diff_b);
    
    pb0 = pb1;
% disp('the total error is')
 diff = diff_v+diff_b;
%  disp(diff);
 if loop>=10000;
     disp(['iteration ',num2str(iter),' not converge']);
     disp(['diff is ',num2str(diff)]);
     break;
 end     
end
disp(['iteration ',num2str(iter),' finished']);
end

polissue = zeros(nzxxb,1);
for i = 1:nzxxb
    tmp = issue(i,:);
    polissue(i) = tmp(polb(i));
end