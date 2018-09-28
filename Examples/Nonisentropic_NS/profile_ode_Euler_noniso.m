function out = profile_ode_Euler_noniso(x,y,s,p)

% nonisentropic ideal polytropic gas

u = y(1);
e = y(2);



out = [ (1/(2*p.mu+p.eta))*(u-1+p.Gamma*(e/u-p.e_minus)); ...
            (1/p.nu)*(e-p.e_minus-(u-1)^2/2+(u-1)*p.Gamma*p.e_minus)];
        
        
        








