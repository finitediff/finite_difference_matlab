function out = fd_jac(U_n,U_o,K,H,p)
%
%Jacobian of finite difference scheme



%Jacobian

mu = p.mu;
nu = p.nu;
u_o = U_o(1,2:end-1);
u_om = U_o(1,1:end-2);
u_op = U_o(1,3:end);
u_n = U_n(1,2:end-1);
u_nm = U_n(1,1:end-2);
u_np = U_n(1,3:end);
v_o = U_o(2,2:end-1);
v_om = U_o(2,1:end-2);
v_op = U_o(2,3:end);
v_n = U_n(2,2:end-1);
v_nm = U_n(2,1:end-2);
v_np = U_n(2,3:end);
%
out{1}{1}{1} = 0;
%
out{1}{1}{2} = (u_n./2 + u_o./2).*(v_n./2 + v_o./2) + 1./K;
%
out{1}{1}{3} = 0;
%
out{1}{2}{1} = 1./(2.*H.^2);
%
out{1}{2}{2} = -mu./2 + (u_n./2 + u_o./2).^2./2 + 3.*(v_n./2 + v_o./2).^2./2 - 1./H.^2;
%
out{1}{2}{3} = 1./(2.*H.^2);
%
out{2}{1}{1} = -1./(2.*H.^2);
%
out{2}{1}{2} = (-u_n./2 - u_o./2).*(u_n./2 + u_o./2) - (u_n./2 + u_o./2).^2./2 - (v_n./2 + v_o./2).^2./2 + 1./2 + H.^(-2);
%
out{2}{1}{3} = -1./(2.*H.^2);
%
out{2}{2}{1} = 0;
%
out{2}{2}{2} = nu + (-u_n./2 - u_o./2).*(v_n./2 + v_o./2) + 1./K;
%
out{2}{2}{3} = 0;
