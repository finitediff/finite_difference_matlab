addpath('../../Core')
clc; clear all; close all; beep off; orig_dir = cd;

% parameters

p.Gamma = 2/3;
p.u_plus = 0.3;
p.mu = 1;
p.eta = -2/3*p.mu;
p.cnu= 1;
p.kappa = 2*p.mu;

% dependent parameters
u_star = p.Gamma/(p.Gamma+2);
alpha = (p.Gamma+2-p.Gamma*p.u_plus)/(p.u_plus-u_star);
p.e_plus = p.u_plus*alpha*(p.u_plus-1)/(p.Gamma*(p.Gamma+2-alpha));
p.e_minus = (p.u_plus-1)*(p.Gamma+2)/(p.Gamma*(p.Gamma+2-alpha));
p.rho_plus = 1/p.u_plus;
p.rho_minus = 1;
p.u_minus = 1;
p.spd = 0;
p.v = 0;
p.nu = p.kappa/p.cnu;
p.R = p.cnu*p.Gamma;
p.k = p.kappa;
p.G = p.Gamma;




%
% profile
%

s.I = 20;
s.R = s.I;
s.L = -s.I;

s.n = 2; % this is the dimension of the profile ode
% we divide the domain in half to deal with the 
% non-uniqueness caused by translational invariance
% s.side = 1 means we are solving the profile on the interval [0,X]
s.side=1; 
s.F=@profile_ode_Euler_noniso; % F is the profile ode
 s.Flinear = @Flinear_noniso;
%s.Flinear = @J; % J is the profile ode Jacobian
s.UL = [p.u_minus;p.e_minus]; % These are the endstates of the profile and its derivative at x = -infty
s.UR = [p.u_plus;p.e_plus]; % These are the endstates of the profile and its derivative at x = +infty
s.phase = 0.5*(s.UL+s.UR); % this is the phase condition for the profile at x = 0
s.order = 1; % this indicates to which componenet the phase conditions is applied
s.stats = 'on'; % this prints data and plots the profile as it is solved
[p,s] = profile_flux(p,s); 

% hold on
% x = linspace(s.L,s.R,200);
% y = zeros(2,length(x));
% for j = 1:length(x)
%     y(:,j)= soln(x(j),s);
% end
% plot(x,y);

dom = linspace(-10,10,300);

bc_L = soln(dom(1),s);
bc_L = [1/bc_L(1);bc_L];
bc_R = soln(dom(end),s);
bc_R = [1/bc_R(1);bc_R];

bc_L_fun = @(U_n,U_o,K,H,p)[U_n(1,1)-bc_L(1);U_n(2,1)-bc_L(2); ...
    U_n(3,1)-bc_L(3)];
bc_R_fun = @(U_n,U_o,K,H,p)[U_n(1,end)-U_n(1,end-1);U_n(2,end)-bc_R(2); ...
    U_n(3,end)-bc_R(3)];
bc_L_jac_fun = @(U_n,U_o,K,H,p)[1 0 0;0 1 0; 0 0 1];
bc_R_jac_fun = @(U_n,U_o,K,H,p)[-1 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];

% -------------------------------------------------------------------------
% Finite Difference Code
% -------------------------------------------------------------------------



U_n = zeros(3,length(dom));
for j = 1:length(dom)
   temp = soln(dom(j),s);
   U_n(1,j) = 1/temp(1);
   U_n(2,j) = temp(1);
   U_n(3,j) = temp(2);
end
U_o = U_n+0.5*sin(dom).*exp(-2*dom.^2);

% delta x
H = dom(2)-dom(1);

% time evolution code tolerance
tol = 1e-8;

figure;
U_orig = U_o;
Ubar = U_n;

times = [0.0001, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0];
tprevious = 0;
for (T = times)
    newT = T - tprevious;
    time_steps = 30;
    K = newT / time_steps;
    tprevious = T;

    p.K = K;
    p.H = H;
    p.T = p.K*time_steps;



    for j = 1:time_steps

        j;

        U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,@fd_F,@fd_jac, ...
        bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);

        %clf;
        %hold on;
        %plot(dom,Ubar,'-r','LineWidth',1);
        %plot(dom,U_orig,'--k','LineWidth',1);
        %plot(dom,U_n,'-g','LineWidth',2);
        %drawnow;

        U_o = U_n;

    end

    close all;
    % 
    % label = {'Plot of \rho','Plot of u','Plot of e'};


    label = {'\rho','u','e'};
    labelClean = {'rho','u','e'};
    file_label = {'rho','u','e'};
    for k = 1:3
        figure
        hold on;
        plot(dom,Ubar(k,:),'-r','LineWidth',1);
        plot(dom,U_orig(k,:),'--k','LineWidth',1);
        plot(dom,U_n(k,:),'-g','LineWidth',2);
        h = xlabel('x');
        set(h,'FontSize',30);
        h = ylabel(label{k});
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
        %fig_name = ['full_gas_T_',num2str(round(p.T)),'_',file_label{k}];
        
        cd('Figures');
        saveas(gca,['full_gas_',labelClean{k},'_',strrep(num2str(T),'.','_'),'.fig']);
        saveas(gca,['full_gas_',labelClean{k},'_',strrep(num2str(T),'.','_'),'.epsc']);
        cd('../');
        

    %     texfig(p,fig_name)

    end
end








