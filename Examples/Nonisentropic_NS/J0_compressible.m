function out = J0_compressible(U,p)

% U is (n X N) where N is the number of nodes and n is the dimension of the
% system.

n = size(U,1);
N = size(U,2);

out = zeros(n,n,N);
out(1,1,:) = ones(N,1);
out(2,1,:) = U(2,:);
out(2,2,:) = U(1,:);
out(3,1,:) = U(3,:)+0.5*U(2,:).^2;
out(3,2,:) = U(1,:).*U(2,:);
out(3,3,:) = U(1,:);