
function jacobian = createJacobian(matrices, time, K, H, parameters)
    V = matrices{1};
    U = matrices{2};
    e = matrices{3};
    cnu = parameters.cnu;
    mu = parameters.mu;
    k = parameters.k;
    Gamma = parameters.Gamma;
    
    %Create the empty matrix.
    jacobianDimensions = length(matrices)*(length(matrices{1}{1})-2);
    jacobian = zeros(jacobianDimensions);
    quadrantDimensions = length(matrices{1}{1})-2;
    n = time -1;
            
    %Loop through quadrant (1,1).
    for row = 1:quadrantDimensions
        for col = 1:quadrantDimensions
            j = col+1;
            if (row - col == 1)
                jacobian(row+quadrantDimensions*0,col+quadrantDimensions*0) = -1/(4*H);
            end
            if (row == col)
                jacobian(row+quadrantDimensions*0,col+quadrantDimensions*0) = 1/K;
            end
            if (row - col == -1)
                jacobian(row+quadrantDimensions*0,col+quadrantDimensions*0) = 1/(4*H);
            end
            
        end
    end
            
            
    %Loop through quadrant (1,2).
    for row = 1:quadrantDimensions
        for col = 1:quadrantDimensions
            j = col+1;
            if (row - col == 1)
                jacobian(row+quadrantDimensions*0,col+quadrantDimensions*1) = 1/(4*H);
            end
            if (row == col)
                jacobian(row+quadrantDimensions*0,col+quadrantDimensions*1) = 0;
            end
            if (row - col == -1)
                jacobian(row+quadrantDimensions*0,col+quadrantDimensions*1) = -1/(4*H);
            end
            
        end
    end
            
            
    %Loop through quadrant (1,3).
    for row = 1:quadrantDimensions
        for col = 1:quadrantDimensions
            j = col+1;
            if (row - col == 1)
                jacobian(row+quadrantDimensions*0,col+quadrantDimensions*2) = 0;
            end
            if (row == col)
                jacobian(row+quadrantDimensions*0,col+quadrantDimensions*2) = 0;
            end
            if (row - col == -1)
                jacobian(row+quadrantDimensions*0,col+quadrantDimensions*2) = 0;
            end
            
        end
    end
            
            
    %Loop through quadrant (2,1).
    for row = 1:quadrantDimensions
        for col = 1:quadrantDimensions
            j = col+1;
            if (row - col == 1)
                jacobian(row+quadrantDimensions*1,col+quadrantDimensions*0) = Gamma*e(n+1,j)/(4*H*V(n+1,j)^2) - mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/(4*H*V(n+1,j)^2);
            end
            if (row == col)
                jacobian(row+quadrantDimensions*1,col+quadrantDimensions*0) = -Gamma*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/V(n+1,j)^2 + 2*Gamma*e(n+1,j)*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^3 + mu*((-2*U(n+1,j) + U(n+1,j+1) + U(n+1,j-1))/(2*H^2) + (-2*U(n,j) + U(n,j+1) + U(n,j-1))/(2*H^2))/V(n+1,j)^2 - 2*mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^3;
            end
            if (row - col == -1)
                jacobian(row+quadrantDimensions*1,col+quadrantDimensions*0) = -Gamma*e(n+1,j)/(4*H*V(n+1,j)^2) + mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/(4*H*V(n+1,j)^2);
            end
            
        end
    end
            
            
    %Loop through quadrant (2,2).
    for row = 1:quadrantDimensions
        for col = 1:quadrantDimensions
            j = col+1;
            if (row - col == 1)
                jacobian(row+quadrantDimensions*1,col+quadrantDimensions*1) = -1/(4*H) - mu*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/(4*H*V(n+1,j)^2) - mu/(2*H^2*V(n+1,j));
            end
            if (row == col)
                jacobian(row+quadrantDimensions*1,col+quadrantDimensions*1) = 1/K + mu/(H^2*V(n+1,j));
            end
            if (row - col == -1)
                jacobian(row+quadrantDimensions*1,col+quadrantDimensions*1) = 1/(4*H) + mu*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/(4*H*V(n+1,j)^2) - mu/(2*H^2*V(n+1,j));
            end
            
        end
    end
            
            
    %Loop through quadrant (2,3).
    for row = 1:quadrantDimensions
        for col = 1:quadrantDimensions
            j = col+1;
            if (row - col == 1)
                jacobian(row+quadrantDimensions*1,col+quadrantDimensions*2) = -Gamma/(4*H*V(n+1,j));
            end
            if (row == col)
                jacobian(row+quadrantDimensions*1,col+quadrantDimensions*2) = -Gamma*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2;
            end
            if (row - col == -1)
                jacobian(row+quadrantDimensions*1,col+quadrantDimensions*2) = Gamma/(4*H*V(n+1,j));
            end
            
        end
    end
            
            
    %Loop through quadrant (3,1).
    for row = 1:quadrantDimensions
        for col = 1:quadrantDimensions
            j = col+1;
            if (row - col == 1)
                jacobian(row+quadrantDimensions*2,col+quadrantDimensions*0) = Gamma*U(n+1,j)*e(n+1,j)/(4*H*V(n+1,j)^2) - U(n+1,j)*mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/(4*H*V(n+1,j)^2) - k*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/(4*H*V(n+1,j)^2*cnu);
            end
            if (row == col)
                jacobian(row+quadrantDimensions*2,col+quadrantDimensions*0) = -Gamma*e(n+1,j)*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/V(n+1,j)^2 + U(n+1,j)*(-Gamma*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/V(n+1,j)^2 + 2*Gamma*e(n+1,j)*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^3) + U(n+1,j)*mu*((-2*U(n+1,j) + U(n+1,j+1) + U(n+1,j-1))/(2*H^2) + (-2*U(n,j) + U(n,j+1) + U(n,j-1))/(2*H^2))/V(n+1,j)^2 - 2*U(n+1,j)*mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^3 + mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))^2/V(n+1,j)^2 + k*((-2*e(n+1,j) + e(n+1,j+1) + e(n+1,j-1))/(2*H^2) + (-2*e(n,j) + e(n,j+1) + e(n,j-1))/(2*H^2))/(V(n+1,j)^2*cnu) - 2*k*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/(V(n+1,j)^3*cnu);
            end
            if (row - col == -1)
                jacobian(row+quadrantDimensions*2,col+quadrantDimensions*0) = -Gamma*U(n+1,j)*e(n+1,j)/(4*H*V(n+1,j)^2) + U(n+1,j)*mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/(4*H*V(n+1,j)^2) + k*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/(4*H*V(n+1,j)^2*cnu);
            end
            
        end
    end
            
            
    %Loop through quadrant (3,2).
    for row = 1:quadrantDimensions
        for col = 1:quadrantDimensions
            j = col+1;
            if (row - col == 1)
                jacobian(row+quadrantDimensions*2,col+quadrantDimensions*1) = -Gamma*e(n+1,j)/(4*H*V(n+1,j)) - U(n+1,j)/(4*H) - U(n+1,j)*mu*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/(4*H*V(n+1,j)^2) + mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/(2*H*V(n+1,j)) - U(n+1,j)*mu/(2*H^2*V(n+1,j));
            end
            if (row == col)
                jacobian(row+quadrantDimensions*2,col+quadrantDimensions*1) = Gamma*((e(n+1,j+1) - e(n+1,j-1))/(4*H) + (e(n,j+1) - e(n,j-1))/(4*H))/V(n+1,j) - Gamma*e(n+1,j)*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2 - mu*((-2*U(n+1,j) + U(n+1,j+1) + U(n+1,j-1))/(2*H^2) + (-2*U(n,j) + U(n,j+1) + U(n,j-1))/(2*H^2))/V(n+1,j) + mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2 + U(n+1,j)/K + (U(n+1,j) - U(n,j))/K + (U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H) + U(n+1,j)*mu/(H^2*V(n+1,j));
            end
            if (row - col == -1)
                jacobian(row+quadrantDimensions*2,col+quadrantDimensions*1) = Gamma*e(n+1,j)/(4*H*V(n+1,j)) + U(n+1,j)/(4*H) + U(n+1,j)*mu*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/(4*H*V(n+1,j)^2) - mu*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/(2*H*V(n+1,j)) - U(n+1,j)*mu/(2*H^2*V(n+1,j));
            end
            
        end
    end
            
            
    %Loop through quadrant (3,3).
    for row = 1:quadrantDimensions
        for col = 1:quadrantDimensions
            j = col+1;
            if (row - col == 1)
                jacobian(row+quadrantDimensions*2,col+quadrantDimensions*2) = -Gamma*U(n+1,j)/(4*H*V(n+1,j)) - 1/(4*H) - k*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/(4*H*V(n+1,j)^2*cnu) - k/(2*H^2*V(n+1,j)*cnu);
            end
            if (row == col)
                jacobian(row+quadrantDimensions*2,col+quadrantDimensions*2) = -Gamma*U(n+1,j)*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/V(n+1,j)^2 + Gamma*((U(n+1,j+1) - U(n+1,j-1))/(4*H) + (U(n,j+1) - U(n,j-1))/(4*H))/V(n+1,j) + 1/K + k/(H^2*V(n+1,j)*cnu);
            end
            if (row - col == -1)
                jacobian(row+quadrantDimensions*2,col+quadrantDimensions*2) = Gamma*U(n+1,j)/(4*H*V(n+1,j)) + 1/(4*H) + k*((V(n+1,j+1) - V(n+1,j-1))/(4*H) + (V(n,j+1) - V(n,j-1))/(4*H))/(4*H*V(n+1,j)^2*cnu) - k/(2*H^2*V(n+1,j)*cnu);
            end
            
        end
    end
            