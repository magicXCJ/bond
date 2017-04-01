nz = Ne*Ne_id*Nj;
nzx = nz*nx;
nzxx = nz*nx*nx;
nzxxb = nb*nz*nx*nx;

egrid = [kron(ones(Ne_id,1),eta),kron(eta_id,ones(Ne,1))];
zgrid =[kron(ones(Nj,1),egrid),kron([0;1],ones(Ne*Ne_id,1))];
zxgrid = [kron(ones(nx,1),zgrid),kron(x,ones(nz,1))];
%z normal,z normal idosyncratic, z jump,x0,x-1
zxxgrid = [kron(ones(nx,1),zxgrid),kron(x,ones(nzx,1))];
%z normal,z normal idosyncratic, z jump,x0,x-1,b0
zxxbgrid = [kron(ones(nb,1),zxxgrid),kron(b,ones(nzxx,1))];

% transition matrix for (x0,x-1) to (x1,x0)
pxx = kron(ones(nx,nx),px).*kron(ones(nx,1),kron(eye(nx),ones(1,nx)));

sig_nzxx = kron(sig,ones(nzx,1));

%cash flow growth (probmatic think again)
%cash flow growth t- to t0
ge0 = exp(mu_e + phi_e*(zxxgrid(:,5) - x_bar)+sig_nzxx.*sig_e.*sig_c.*zxxgrid(:,1)-0.5*sig_nzxx.^2*sig_e^2*sig_c^2  ...
      +sig_nzxx.*sig_eid.*zxxgrid(:,2)-0.5*sig_nzxx.^2*sig_eid^2 ...
       +(d_e.*zxxgrid(:,2) - (exp(d_e)-1)*lambda_e*zxxgrid(:,5)));
   
% %cash flow growth t0 to t1

% prob to e,e_id;
w_e = kron(w_id,w);
%prob of transition to z condition on x-1   
pz = kron(ones(Nj*nx*nx,1),w_e).*(1-lambda_e*zxxgrid(:,5)-zxxgrid(:,3)+2*zxxgrid(:,3).*lambda_e.*zxxgrid(:,5));

% transition matrix for zxx;
pzxx = kron(pxx,ones(nz,nz)).* kron(ones(nzxx,1),pz');

% construct pricing kernel
Mxx = kron(ones(nx,nx),M).*kron(ones(nx,1),kron(eye(nx),ones(1,nx)));
Mz = exp(-gamma.*sig_nzxx.*sig_c.*zxxgrid(:,1));
Mzxx = kron(Mxx,ones(nz,nz)).* kron(ones(nzxx,1),Mz');

% risk neutural prob;
pzxx_q = pzxx.*Mzxx;
% risk neutural form x0 to zxx;
px_zxx_q = zeros(nx,nzxx);
for i = 1:nx
    px_zxx_q(i,:) = pzxx_q((i-1)*nz+1,:);
end

% store expanded pe
pe_zxx = kron(ones(nx,1),kron(pe,ones(nz,1)));

%bond issurance;
issue = kron(ones(nzxxb,1),b')-(1-kappa)*kron(ones(1,nb),zxxbgrid(:,6)./kron(ones(nb,1),ge0));
ne_penalty = (issue<0).*issue*10000000;
%inertia : the smallest positive bond issurance, fixed issurance cost = 0
[tmp] = min(abs(issue).*(issue>0),[],2);
[~,inertia] = min(abs(issue - kron(ones(1,nb),tmp)),[],2);
adjcost = phi_0*(issue>kron(ones(1,nb),tmp))+phi_1*issue;
div_cons = (1-tau) - ((1-tau)*coupon+kappa)*kron(ones(1,nb),zxxbgrid(:,6)./kron(ones(nb,1),ge0))-adjcost;                
div_cons = div_cons+ne_penalty;

%store default payoff
defpayoff = kron(ones(1,nb),ge0).*(1-alpha).*kron(ones(1,nb),pe_zxx)./kron(ones(nzxx,1),b');
%defpayoff = min(1+coupon,defpayoff);

% shock, state t, state t-1, b_t/e_(t-1)
pv = zeros(nzxxb,1);
%pv = zeros(nz,nx,nx,nb);
pb = ones(nx,nb)*1;
% 
pb0 = pb;
pb1 = pb0;
pb0=pb;

pv0 = pv;
pv1 = pv0;
polb = pv0;



%next period bond level

diff_b = 1;
while diff_b> sqrt(eps)%0.0001
diff_v = 1;
while diff_v > sqrt(eps)%0.0001    
    div = div_cons + kron(ones(nx*nb,1),kron(pb0,ones(nz,1))).*issue;
    expv = pzxx_q * (kron(ge0,ones(1,nb)).*reshape(pv0,nzxx,nb));
    [vmax, polb] = max(div+kron(ones(nb,1),expv),[],2);
    %pv1 = vmax;
    pv1 = max(0,vmax); 


    diff_v = norm(pv1(:) - pv0(:));
    %diff_v = sumabs(pv1(:) - pv0(:));
    %disp(diff_v);
    pv0 = pv1;
end
polbreal = b(polb);
% given default policy, calculate the bond price, policy consistent
%disp(pv0)

pv0_r = reshape(pv0,nzxx,nb);

Id = (pv0_r<=0);
ndef = (1-Id)*(coupon+1);
def = Id.*defpayoff;
pb1 = px_zxx_q*(ndef+def);
% pb_zxx = pzxx_q*(ndef+def);
% pb_r = reshape(pb_zxx,nz,nx,nx,nb);
% 
% pb1 = squeeze(pb_r(1,:,1,:));
diff_b = norm(pb1(:) - pb0(:));
disp('the bond function error is')
disp(diff_b);
    
pb0 = pb1;
end

pb_kappa1 = pb0;
pv_kappa1 = pv0;
polb_kappa1 = polb;
polbreal_kapp1 = polbreal;