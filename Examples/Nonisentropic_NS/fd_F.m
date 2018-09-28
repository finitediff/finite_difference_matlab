function out = fd_F(U_n,U_o,K,H,p)
%
%Function of the finite difference scheme

nu = p.nu;
mu = p.mu;
eta = p.eta;
cnu = p.cnu;
G = p.G;
rho_o = U_o(1,2:end-1);
rho_om = U_o(1,1:end-2);
rho_op = U_o(1,3:end);
rho_n = U_n(1,2:end-1);
rho_nm = U_n(1,1:end-2);
rho_np = U_n(1,3:end);
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


out = [rho_n.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H)) + u_n.*((-rho_nm + rho_np)./(4.*H) + (-rho_om + rho_op)./(4.*H)) + (rho_n - rho_o)./K;
G.*rho_n.*((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H)) + 2.*rho_n.*u_n.*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H)) - (eta + 2.*mu).*((-2.*u_n + u_nm + u_np)./(2.*H.^2) + (-2.*u_o + u_om + u_op)./(2.*H.^2)) + (G.*e_n + u_n.^2).*((-rho_nm + rho_np)./(4.*H) + (-rho_om + rho_op)./(4.*H)) + rho_n.*(u_n - u_o)./K + u_n.*(rho_n - rho_o)./K;
-nu.*((-2.*e_n + e_nm + e_np)./(2.*H.^2) + (-2.*e_o + e_om + e_op)./(2.*H.^2)) - u_n.*(eta + 2.*mu).*((-2.*u_n + u_nm + u_np)./(2.*H.^2) + (-2.*u_o + u_om + u_op)./(2.*H.^2)) - (eta + 2.*mu).*((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H)).^2 + ((-e_nm + e_np)./(4.*H) + (-e_om + e_op)./(4.*H)).*(G.*rho_n.*u_n + rho_n.*u_n) + ((-rho_nm + rho_np)./(4.*H) + (-rho_om + rho_op)./(4.*H)).*(G.*e_n.*u_n + u_n.*(e_n + u_n.^2./2)) + ((-u_nm + u_np)./(4.*H) + (-u_om + u_op)./(4.*H)).*(G.*e_n.*rho_n + rho_n.*u_n.^2 + rho_n.*(e_n + u_n.^2./2)) + rho_n.*u_n.*(u_n - u_o)./K + rho_n.*(e_n - e_o)./K + (e_n + u_n.^2./2).*(rho_n - rho_o)./K;
];

