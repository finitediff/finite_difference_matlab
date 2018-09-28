addpath('../../Core')
clc; clear all; close all; beep off;

num_inf = 30;

H_0 = 0.1;
K_0 = 0.1;


T_0 = 3; %T = 0, T = 0.1, T = 0.5, T = 1, T = 10, T = 3


eps_vals = 1;

for cnt = 1:length(eps_vals)


    eps = eps_vals(cnt);
    
    %Set the domain variables.
    K = eps^2*K_0;
    H = eps*H_0;
    s.L = -num_inf;
    s.R = num_inf;
    dom = linspace(s.L,s.R,1+round((2*num_inf)/H));
    H = (dom(2)-dom(1));

    %Write the boundary functions
    bc_L_fun = @(U_n,U_o,K,H,p)[U_n(1)-2];
    bc_R_fun = @(U_n,U_o,K,H,p)[U_n(end)];
    bc_L_jac_fun = @(U_n,U_o,K,H,p)(1);
    bc_R_jac_fun = @(U_n,U_o,K,H,p)(1);

    % -------------------------------------------------------------------------
    % Finite Difference Code
    % -------------------------------------------------------------------------

    % time steps
    time_steps = round(T_0/K)

    % delta x and delta y

    startFun = @(x)(1-tanh(x./2));
    U_o = startFun(dom);
    U_orig = U_o;
    U_o = U_o+1*sin(dom).*exp(-0.2*(dom+0).^2);
    U_n = U_o;

    % time evolution code tolerance
    tol = 1e-8;

    figure;

    Ubar = U_n;

    p = {};

    for j = 1:time_steps

       
        U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,@fd_F,@fd_jac, ...
        bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);

        clf;
        hold on;
        plot(dom,Ubar,'-r','LineWidth',1);
        plot(dom,U_orig,'--k','LineWidth',1);
        plot(dom,U_n,'-g','LineWidth',2);
        axis([s.L,s.R,-.2,2.4]);
        h = xlabel('x');
        set(h,'FontSize',18);
        h = gca;
        set(h,'FontSize',18);
        h = ylabel('u');
        set(h,'FontSize',18);
        drawnow;


        U_o = U_n;

    end

end







    
