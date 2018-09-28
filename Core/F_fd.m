function out = F_fd(U,old_U,fun_J0,fun_J1,fun_B,fun_Bx,p,delt,delx)

n = size(U,1);
N = size(U,2);

J0 = fun_J0(U,p);
J1 = fun_J1(U,p);
BU = fun_B(U,p);
BUx = fun_Bx(U,p);


U_t = (U(:,2:end-1)-old_U(:,2:end-1))/delt;
U_xx = 0.5*(U(:,1:end-2)-2*U(:,2:end-1)+U(:,3:end))/delx^2 + ...
        0.5*(old_U(:,1:end-2)-2*old_U(:,2:end-1)+old_U(:,3:end))/delx^2;
U_x = 0.5*(U(:,1:end-2)-U(:,3:end))/(2*delx) + ...
        0.5*(old_U(:,1:end-2)-old_U(:,3:end))/(2*delx);

    

    
out = mat_mult_vec(J0(:,:,2:end-1),U_t) + ...
    mat_mult_vec(J1(:,:,2:end-1),U_x) + ...
    mat_mult_vec(BU(:,:,2:end-1),U_xx)+ ...
    mat_mult_vec(BUx(:,:,2:end-1),U_x);


