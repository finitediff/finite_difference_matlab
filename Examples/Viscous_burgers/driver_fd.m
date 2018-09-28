addpath('../../Core')
clc; clear all; close all; beep off;

% mathcalE = 2.7; % stable wave
mathcalE = 7; %


%Set the domain variables.
s.L = -20;
s.R = 20;
dom = linspace(s.L,s.R,1600);


%Write the boundary functions
bc_L_fun = @(U_n,U_o,K,H,p)[U_n(1)-2];
bc_R_fun = @(U_n,U_o,K,H,p)[U_n(end)];
bc_L_jac_fun = @(U_n,U_o,K,H,p)(1);
bc_R_jac_fun = @(U_n,U_o,K,H,p)(1);

% -------------------------------------------------------------------------
% Finite Difference Code
% -------------------------------------------------------------------------

%U_o = U_n+0.03*exp(-2*(dom+2).^2);

% time steps
time_steps = 1000;


% delta x and delta y
H = (dom(2)-dom(1))
K = (.1)

startFun = @(x)[1-tanh(x./2)]
U_orig = startFun(dom)
U_o = U_orig+4.5.*sin(.3*dom).*exp(-2*(dom).^2);
U_n = U_o;

% time evolution code tolerance
tol = 1e-8;
p = {}

figure;
Ubar = U_n;

for j = 1:time_steps
    
    j
    
    U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,@fd_F,@fd_jac, ...
    bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);

    clf;
    hold on;
    plot(dom,Ubar,'-r','LineWidth',1);
    plot(dom,U_orig,'--k','LineWidth',1);
    plot(dom,U_n,'-g','LineWidth',2);
    axis([.3*s.L,.3*s.R,-.1,2.1]);
    drawnow;


    U_o = U_n;
    
end















    
