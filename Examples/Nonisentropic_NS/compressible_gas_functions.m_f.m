
%loop through equations, creating a vector for all j's
function outVector = f(matrices, parameters, K, H, n)
    V = matrices{1};
    U = matrices{2};
    e = matrices{3};
    cnu = parameters.cnu;
    mu = parameters.mu;
    k = parameters.k;
    Gamma = parameters.Gamma;
    ind = length(matrices{0}{0})-1;
    outVector = zeros(3*ind,1);
    
    %Add values for equation 1
        for j = 0*ind+1:1*ind
            outVector(j) = (V(n+1,j) - V(n,j))/K - (U(n+1,j+1) - U(n+1,j-1))/(4*H) - (U(n,j+1) - U(n,j-1))/(4*H) + (V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H);
        end
        
    %Add values for equation 2
        for j = 1*ind+1:2*ind
            outVector(j) = Gamma*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/V(n+1,j) - Gamma*e(n+1,j)*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2 - mu*((-2*U(n+1,j) + U(n+1,j+1) + U(n+1,j-1))/(2*H^2) + (-2*U(n,j) + U(n,j+1) + U(n,j-1))/(2*H^2))/V(n+1,j) + mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2 + (U(n+1,j) - U(n,j))/K + (U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H);
        end
        
    %Add values for equation 3
        for j = 2*ind+1:3*ind
            outVector(j) = Gamma*e(n+1,j)*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/V(n+1,j) + U(n+1,j)*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H)) + U(n+1,j)*(Gamma*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/V(n+1,j) - Gamma*e(n+1,j)*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2) - U(n+1,j)*mu*((-2*U(n+1,j) + U(n+1,j+1) + U(n+1,j-1))/(2*H^2) + (-2*U(n,j) + U(n,j+1) + U(n,j-1))/(2*H^2))/V(n+1,j) + U(n+1,j)*mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2 - mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))^2/V(n+1,j) - k*((-2*e(n+1,j) + e(n+1,j+1) + e(n+1,j-1))/(2*H^2) + (-2*e(n,j) + e(n,j+1) + e(n,j-1))/(2*H^2))/(V(n+1,j)*cnu) + k*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/(V(n+1,j)^2*cnu) + U(n+1,j)*(U(n+1,j) - U(n,j))/K + (e(n+1,j) - e(n,j))/K + (e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H);
        end
        
    

    