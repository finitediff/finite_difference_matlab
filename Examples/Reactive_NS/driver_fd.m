addpath('../../Core')
clc; clear all; close all; beep off;

orig_dir = cd;
folder = 'rNS';
subfolder = 'data';

% parameters

epsilon = 0.1;
Ti_weight = 0.99;
num_inf = 30;

%mathcalE = 2.7; % stable wave
mathcalE = 7; % unstable wave, but not every perturbation will show that quickly

% file_name = ['rts_E_',num2str(10000000000*mathcalE),'eps',num2str(10000000000*epsilon)];
% ld = retrieve_it(orig_dir,folder,file_name,subfolder);
file_name = ['profile_E_',num2str(10000000000*mathcalE),'eps',num2str(10000000000*epsilon)];
ld = load(['data/',file_name]);
d = ld.var;
p = d{1};

% return

% ld = retrieve_it(orig_dir,'rNS',file_name,'data');
% p = ld.var{1};
s = ld.var{2};
p.cnu = 1;
dom = linspace(-num_inf,num_inf,1600);

temp = soln(s.L,s);
bc_L = [1-temp(1);temp(1);temp(2);temp(3)];

% bc_L = [1/bc_L(1);bc_L];
temp = soln(s.R,s);
bc_R = [1-temp(1);temp(1);temp(2);temp(3)];

bc_L_fun = @(U_n,U_o,K,H,p)[U_n(1,1)-bc_L(1);U_n(2,1)-bc_L(2); ...
    U_n(3,1)-bc_L(3);U_n(4,1)-bc_L(4)];
bc_R_fun = @(U_n,U_o,K,H,p)[U_n(1,end)-bc_R(1);U_n(2,end)-bc_R(2); ...
    U_n(3,end)-bc_R(3);U_n(4,end)-U_n(4,end-1)];
bc_L_jac_fun = @(U_n,U_o,K,H,p)[1,0,0,0; 0,1,0,0; 0,0,1,0; 0,0,0,1];
bc_R_jac_fun = @(U_n,U_o,K,H,p)[0 0 0 0 1,0,0,0; 0 0 0 0 0,1,0,0;...
 0 0 0 0 0,0,1,0; 0 0 0 -1 0,0,0,1];

% -------------------------------------------------------------------------
% Finite Difference Code
% -------------------------------------------------------------------------



U_n = zeros(4,length(dom));
for j = 1:length(dom)
    if dom(j) >= s.L
        if dom(j) <= s.R
            temp = soln(dom(j),s);
        else
            temp = soln(s.R,s);
        end
    else
        temp = soln(s.L,s);
    end
   U_n(:,j) = [1-temp(1);temp(1);temp(2);temp(3)];
end

% return
U_o = U_n+0.015*sin(dom).*exp(-4*(dom+2).^2);

% delt T
K = 0.03;

% time steps
time_steps = 100;

% delta x
H = dom(2)-dom(1)

% time evolution code tolerance
tol = 1e-8;

figure;
U_orig = U_o;
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
    axis([-.5*num_inf,0.1*num_inf,-0.05,1.05]);
    drawnow;


    U_o = U_n;

end











    
