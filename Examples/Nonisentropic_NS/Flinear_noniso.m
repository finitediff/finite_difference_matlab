function out = Flinear_noniso(w,p)

u = w(1);
e = w(2);

out = [(1/(2*p.mu+p.eta))*(1-p.Gamma*(e/u^2)), (1/(2*p.mu+p.eta))*p.Gamma/u;
    -(1/p.nu)*(u-1)+p.Gamma*p.e_minus, (1/p.nu)];


