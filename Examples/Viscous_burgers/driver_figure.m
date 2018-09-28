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

% delta x
H = (dom(2)-dom(1))

% set the initial functions.
startFun = @(x)[1-tanh(x./2)]
U_orig = startFun(dom)
U_o = U_orig+4.5.*sin(.3*dom).*exp(-2*(dom).^2);
U_n = U_o;
Ubar = U_n;

% time evolution code tolerance
tol = 1e-8;
p = {}
figure;

%Begin the time loop
times = [0.0001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 5.0];
tprevious = 0;
for (T = times)
    newT = T - tprevious;
    time_steps = 30;
    K = newT / time_steps;
    tprevious = T;

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

    k = 1;
    label = 'u';
    figure
    hold on;
    plot(dom,Ubar(k,:),'-r','LineWidth',1);
    plot(dom,U_orig(k,:),'--k','LineWidth',1);
    plot(dom,U_n(k,:),'-g','LineWidth',2);
    axis([.3*s.L,.3*s.R,-.1,2.1]);
    h = xlabel('x');
    set(h,'FontSize',30);
    h = ylabel(label);
    %     h = title(label{k});
    %     h = legend(label{k},'Location','Best');
    set(h,'FontSize',30);
    h = gca;
    set(h,'FontSize',24);
    scale = 0.1;
    pos = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-1.5*scale)*pos(4);
    pos(1) = pos(1)+scale*pos(3);
    pos(3) = (1-scale)*pos(3);
    set(gca, 'Position', pos)
    drawnow;
    
    cd('Figures');
    saveas(gca,['burgers_u_',strrep(num2str(T),'.','_'),'.fig']);
    saveas(gca,['burgers_u_',strrep(num2str(T),'.','_'),'.epsc']);
    cd('../');
end










    
