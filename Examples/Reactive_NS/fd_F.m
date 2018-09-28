function out = fd_F(U_n,U_o,K,H,p)
%
%Function of the finite difference scheme

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
for j = 1:length(e_n)
   if cnu*e_n(j)-Ti > 0
       exp_A(j) = exp(-A./(-Ti+e_n(j)/cnu));
   else
       exp_A(j) = 0;
   end
end



out = [-spd.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H)) + (tau_n - tau_o)./K - (-u_nm + u_np)./(4.*H) - (-u_om + u_op)./(4.*H);
-G.*e_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^2 + G.*((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H))./tau_n - nu.*((-2.*u_n + u_nm + u_np)./(2.*H.^2) + (-2.*u_o + u_om + u_op)./(2.*H.^2))./tau_n + nu.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H)).*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./tau_n.^2 - spd.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H)) + (u_n - u_o)./K;
-G.*e_n.*u_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^2 - k.*q.*z_n.*exp_A - nu.*u_n.*((-2.*u_n + u_nm + u_np)./(2.*H.^2) + (-2.*u_o + u_om + u_op)./(2.*H.^2))./tau_n + ((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H)).*(G.*u_n./tau_n - spd) + ((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H)).*(G.*e_n./tau_n - spd.*u_n) - ((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H)).*(nu.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H))./tau_n - nu.*u_n.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./tau_n.^2) - kappa.*((-2.*e_n + e_nm + e_np)./(2.*H.^2) + (-2.*e_o + e_om + e_op)./(2.*H.^2))./(cnu.*tau_n) + kappa.*((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H)).*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H))./(cnu.*tau_n.^2) + u_n.*(u_n - u_o)./K + (e_n - e_o)./K;
-D.*((-2.*z_n + z_nm + z_np)./(2.*H.^2) + (-2.*z_o + z_om + z_op)./(2.*H.^2))./tau_n.^2 + 2.*D.*((-tau_nm + tau_np)./(4.*H) + (-tau_om + tau_op)./(4.*H)).*((-z_nm + z_np)./(4.*H) + (-z_om + z_op)./(4.*H))./tau_n.^3 + k.*z_n.*exp_A - spd.*((-z_nm + z_np)./(4.*H) + (-z_om + z_op)./(4.*H)) + (z_n - z_o)./K;
];

