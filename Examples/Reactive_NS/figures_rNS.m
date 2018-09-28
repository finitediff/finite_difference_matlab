addpath('../../Core')
clc; clear all; close all; beep off;

%   NEEDS SOME WORK, BUT THIS SEEMS TO CONVERGE SUFFICIENTLY TO A TRANSLATE

num_inf = 30;
num_solved = 8.9566;

H_0 = 0.033608236638813;
K_0 = 0.1;
T_0 = 12;
% eps_vals = [1,1/2,1/4,1/8];
eps_vals = 1;%[1,1/2]
for cnt = 1:length(eps_vals)


    eps = eps_vals(cnt);
    
    %Set the domain variables.
    K = eps^2*K_0;
    H = eps*H_0;
    s.L = -num_solved;
    s.R = num_solved;
    %dom = linspace(-8,8,1+round((2*num_inf)/H));

  
    %H = (dom(2)-dom(1));
    orig_dir = cd;
    folder = 'rNS';
    subfolder = 'data';

    % parameters

    epsilon = 0.1;
    Ti_weight = 0.99;

    mathcalE = 2.7; % stable wave
%     mathcalE = 7; % unstable wave, but not every perturbation will show that quickly

    % file_name = ['rts_E_',num2str(10000000000*mathcalE),'eps',num2str(10000000000*epsilon)];
    % ld = retrieve_it(orig_dir,folder,file_name,subfolder);
    file_name = ['profile_E_',num2str(10000000000*mathcalE),'eps',num2str(10000000000*epsilon)];
    ld = load("data/"+file_name);
    d = ld.var;
    p = d{1};
    p.cnu = 1;
    s = ld.var{2};
    
    num_inf = s.R;
    num_inf_L = -2*num_inf;
    num_inf_R = num_inf;
    dom = linspace(num_inf_L,num_inf_R,1+round((num_inf_R-num_inf_L)/H));

    % ld = retrieve_it(orig_dir,'rNS',file_name,'data');
    % p = ld.var{1};
    
    


    %Set up the boundary conditions.
    temp = soln(s.L,s);
    bc_L = [1-temp(1);temp(1);temp(2);temp(3)];
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
    
    % Extend the end of the profile.
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

    hold on;
    plot(dom,U_n,'-k');
    plot(dom,U_n(1,:),'-g');
    % plot(dom,p.Ti*ones(size(U_n(3,:))),'-r');

    % Add perturbation.
    U_o = U_n+0.015*sin(dom).*exp(-4*(dom+2).^2);

    % time steps
    time_steps = round(T_0/K);

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
        plot(dom,Ubar,'-r','LineWidth',2);
        plot(dom,U_orig,'--k','LineWidth',2);
        plot(dom,U_n,'-g','LineWidth',2);
        axis([-3*num_inf,0.5*num_inf,-0.05,1.05]);
        drawnow;


        U_o = U_n;

    end
    waves{cnt}.dom = dom;
    waves{cnt}.U = U_n;

end
n = length(eps_vals);

for j = 1:n
   
    scl = 2^(n-j);
    max_diff = 0;
    
    for k = 1:4
        max_diff = max(max_diff,max(abs(waves{end}.U(k,1:scl:end)-waves{j}.U(k,:))));
    end
    max_diff
    
end

    
