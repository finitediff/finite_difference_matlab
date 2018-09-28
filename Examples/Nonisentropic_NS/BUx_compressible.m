function out = BUx_compressible(U,p)

% U is (n X N) where N is the number of nodes and n is the dimension of the
% system.

n = size(U,1);
N = size(U,2);

out = zeros(n,n,N);

out(3,2,:) = (2*p.mu+p.eta)*ones(N,1);
