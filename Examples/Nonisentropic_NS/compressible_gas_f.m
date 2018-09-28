
%loop through equations, creating a vector for all j's
function outVector = compressible_gas_f(matrices, parameters, K, H, n)
    V = matrices{1};
    U = matrices{2};
    e = matrices{3};
    cnu = parameters.cnu;
    Gamma = parameters.Gamma;
    mu = parameters.mu;
    k = parameters.k;
    ind = length(matrices{1})-1;
    outVector = zeros(3*(ind-2),1);
    
    index_count = 1;
    
    %Add values for equation 1
        for j = 2:ind-1
            outVector(index_count) = (V(n+1,j) - V(n,j))/K - (U(n+1,j+1) - U(n+1,j-1))/(4*H) - (U(n,j+1) - U(n,j-1))/(4*H) + (V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H);
            index_count = index_count + 1;
        end
        
    %Add values for equation 2
        for j = 2:ind-1
            outVector(index_count) = Gamma*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/V(n+1,j) - Gamma*e(n+1,j)*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2 - mu*((-2*U(n+1,j) + U(n+1,j+1) + U(n+1,j-1))/(2*H^2) + (-2*U(n,j) + U(n,j+1) + U(n,j-1))/(2*H^2))/V(n+1,j) + mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2 + (U(n+1,j) - U(n,j))/K + (U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H);
            index_count = index_count + 1;
        end
        
    %Add values for equation 3
        for j = 2:ind-1
            outVector(index_count) = Gamma*e(n+1,j)*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/V(n+1,j) + U(n+1,j)*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H)) + U(n+1,j)*(Gamma*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/V(n+1,j) - Gamma*e(n+1,j)*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2) - U(n+1,j)*mu*((-2*U(n+1,j) + U(n+1,j+1) + U(n+1,j-1))/(2*H^2) + (-2*U(n,j) + U(n,j+1) + U(n,j-1))/(2*H^2))/V(n+1,j) + U(n+1,j)*mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2 - mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))^2/V(n+1,j) - k*((-2*e(n+1,j) + e(n+1,j+1) + e(n+1,j-1))/(2*H^2) + (-2*e(n,j) + e(n,j+1) + e(n,j-1))/(2*H^2))/(V(n+1,j)*cnu) + k*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/(V(n+1,j)^2*cnu) + U(n+1,j)*(U(n+1,j) - U(n,j))/K + (e(n+1,j) - e(n,j))/K + (e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H);
            index_count = index_count + 1;
        end
        
    

    