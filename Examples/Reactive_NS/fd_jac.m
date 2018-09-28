function out = fd_jac(U_n,U_o,K,H,p)
%
%Jacobian of finite difference scheme



%Jacobian

nu = p.nu;
kappa = p.kappa;
D = p.d;
cnu = p.cnu;
q = p.q;
k = p.k;
G = p.g;
A = p.mathcalE;
Ti = p.Ti;
spd = p.s;




tau_o = U_o(1,2:end-1);
tau_om = U_o(1,1:end-2);
tau_op = U_o(1,3:end);
tau_n = U_n(1,2:end-1);
tau_nm = U_n(1,1:end-2);
tau_np = U_n(1,3:end);
u_o = U_o(2,2:end-1);
u_om = U_o(2,1:end-2);
u_op = U_o(2,3:end);
u_n = U_n(2,2:end-1);
u_nm = U_n(2,1:end-2);
u_np = U_n(2,3:end);
e_o = U_o(3,2:end-1);
e_om = U_o(3,1:end-2);
e_op = U_o(3,3:end);
e_n = U_n(3,2:end-1);
e_nm = U_n(3,1:end-2);
e_np = U_n(3,3:end);
z_o = U_o(4,2:end-1);
z_om = U_o(4,1:end-2);
z_op = U_o(4,3:end);
z_n = U_n(4,2:end-1);
z_nm = U_n(4,1:end-2);
z_np = U_n(4,3:end);

expA = zeros(size(e_n));
d_expA = zeros(size(e_n));
for j = 1:length(e_n)
   if cnu*e_n(j)-Ti > 0
       exp_A(j) = exp(-A./(-Ti+e_n(j)/cnu)); 
%        d_expA(j) = expA(j)./(-Ti + cnu.*e_n(j)).^2;% KEEP AS IS 
   else
       exp_A(j) = 0;
%        d_expA(j) = 0;
   end
end



%
out{1}{1}{1} = spd./(4.*H);
%
out{1}{1}{2} = 1./K;
%
out{1}{1}{3} = -spd./(4.*H);
%
out{1}{2}{1} = 1./(4.*H);
%
out{1}{2}{2} = 0;
%
out{1}{2}{3} = -1./(4.*H);
%
out{1}{3}{1} = 0;
%
out{1}{3}{2} = 0;
%
out{1}{3}{3} = 0;
%
out{1}{4}{1} = 0;
%
out{1}{4}{2} = 0;
%
out{1}{4}{3} = 0;
%
out{2}{1}{1} = G.*e_n./(4.*H.*tau_n.^2) - nu.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./(4.*H.*tau_n.^2);
%
out{2}{1}{2} = 2.*G.*e_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^3 - G.*((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H))./tau_n.^2 + nu.*((-2.*u_n + u_nm + u_np)./(2.*H.^2) + (-2.*u_o + u_om + u_op)./(2.*H.^2))./tau_n.^2 - 2.*nu.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H)).*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./tau_n.^3;
%
out{2}{1}{3} = -G.*e_n./(4.*H.*tau_n.^2) + nu.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./(4.*H.*tau_n.^2);
%
out{2}{2}{1} = -nu.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./(4.*H.*tau_n.^2) + spd./(4.*H) - nu./(2.*H.^2.*tau_n);
%
out{2}{2}{2} = 1./K + nu./(H.^2.*tau_n);
%
out{2}{2}{3} = nu.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./(4.*H.*tau_n.^2) - spd./(4.*H) - nu./(2.*H.^2.*tau_n);
%
out{2}{3}{1} = -G./(4.*H.*tau_n);
%
out{2}{3}{2} = -G.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^2;
%
out{2}{3}{3} = G./(4.*H.*tau_n);
%
out{2}{4}{1} = 0;
%
out{2}{4}{2} = 0;
%
out{2}{4}{3} = 0;
%
out{3}{1}{1} = G.*e_n.*u_n./(4.*H.*tau_n.^2) - nu.*u_n.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./(4.*H.*tau_n.^2) - kappa.*((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H))./(4.*H.*cnu.*tau_n.^2);
%
out{3}{1}{2} = -G.*e_n.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./tau_n.^2 + 2.*G.*e_n.*u_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^3 - G.*u_n.*((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H))./tau_n.^2 + nu.*u_n.*((-2.*u_n + u_nm + u_np)./(2.*H.^2) + (-2.*u_o + u_om + u_op)./(2.*H.^2))./tau_n.^2 - ((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H)).*(-nu.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./tau_n.^2 + 2.*nu.*u_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^3) + kappa.*((-2.*e_n + e_nm + e_np)./(2.*H.^2) + (-2.*e_o + e_om + e_op)./(2.*H.^2))./(cnu.*tau_n.^2) - 2.*kappa.*((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H)).*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./(cnu.*tau_n.^3);
%
out{3}{1}{3} = -G.*e_n.*u_n./(4.*H.*tau_n.^2) + nu.*u_n.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./(4.*H.*tau_n.^2) + kappa.*((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H))./(4.*H.*cnu.*tau_n.^2);
%
out{3}{2}{1} = nu.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./(4.*H.*tau_n) - (G.*e_n./tau_n - spd.*u_n)./(4.*H) + (nu.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./tau_n - nu.*u_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^2)./(4.*H) - nu.*u_n./(2.*H.^2.*tau_n);
%
out{3}{2}{2} = -G.*e_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^2 + G.*((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H))./tau_n - nu.*((-2.*u_n + u_nm + u_np)./(2.*H.^2) + (-2.*u_o + u_om + u_op)./(2.*H.^2))./tau_n + nu.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H)).*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./tau_n.^2 - spd.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H)) + u_n./K + (u_n - u_o)./K + nu.*u_n./(H.^2.*tau_n);
%
out{3}{2}{3} = -nu.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./(4.*H.*tau_n) + (G.*e_n./tau_n - spd.*u_n)./(4.*H) - (nu.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./tau_n - nu.*u_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^2)./(4.*H) - nu.*u_n./(2.*H.^2.*tau_n);
%
out{3}{3}{1} = -(G.*u_n./tau_n - spd)./(4.*H) - kappa.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./(4.*H.*cnu.*tau_n.^2) - kappa./(2.*H.^2.*cnu.*tau_n);
%
out{3}{3}{2} = -A.*k.*q.*z_n.*exp_A./(cnu.*(-Ti + e_n./cnu).^2) + G.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./tau_n - G.*u_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^2 + 1./K + kappa./(H.^2.*cnu.*tau_n);
%
out{3}{3}{3} = (G.*u_n./tau_n - spd)./(4.*H) + kappa.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./(4.*H.*cnu.*tau_n.^2) - kappa./(2.*H.^2.*cnu.*tau_n);
%
out{3}{4}{1} = 0;
%
out{3}{4}{2} = -k.*q.*exp_A;
%
out{3}{4}{3} = 0;
%
out{4}{1}{1} = -D.*((-z_nm + z_np)./(4.*H) + (-z_om + z_op)./(4.*H))./(2.*H.*tau_n.^3);
%
out{4}{1}{2} = 2.*D.*((-2.*z_n + z_nm + z_np)./(2.*H.^2) + (-2.*z_o + z_om + z_op)./(2.*H.^2))./tau_n.^3 - 6.*D.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H)).*((-z_nm + z_np)./(4.*H) + (-z_om + z_op)./(4.*H))./tau_n.^4;
%
out{4}{1}{3} = D.*((-z_nm + z_np)./(4.*H) + (-z_om + z_op)./(4.*H))./(2.*H.*tau_n.^3);
%
out{4}{2}{1} = 0;
%
out{4}{2}{2} = 0;
%
out{4}{2}{3} = 0;
%
out{4}{3}{1} = 0;
%
out{4}{3}{2} = A.*k.*z_n.*exp_A./(cnu.*(-Ti + e_n./cnu).^2);
%
out{4}{3}{3} = 0;
%
out{4}{4}{1} = -D.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./(2.*H.*tau_n.^3) - D./(2.*H.^2.*tau_n.^2) + spd./(4.*H);
%
out{4}{4}{2} = D./(H.^2.*tau_n.^2) + k.*exp_A + 1./K;
%
out{4}{4}{3} = D.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./(2.*H.*tau_n.^3) - D./(2.*H.^2.*tau_n.^2) - spd./(4.*H);
