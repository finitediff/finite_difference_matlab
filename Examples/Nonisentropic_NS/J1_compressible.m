function out = J1_compressible(U,p)

% U is (n X N) where N is the number of nodes and n is the dimension of the
% system.

n = size(U,1);
N = size(U,2);

out = zeros(n,n,N);

out(1,1,:) = U(2,:);
out(1,2,:) = U(1,:);
out(2,1,:) =U(2,:).^2+p.Gamma*U(3,:);
out(2,2,:) = 2*U(1,:).*U(2,:);
out(2,3,:) = p.Gamma*U(1,:);
out(3,1,:) = U(2,:).*(U(3,:)+0.5*U(2,:).^2)+p.Gamma*U(2,:).*U(3,:);
out(3,2,:) = U(1,:).*(U(3,:)+0.5*U(2,:).^2)+U(1,:).*U(2,:).^2+p.Gamma*U(1,:).*U(3,:);
out(3,3,:) = U(1,:).*U(2,:)+p.Gamma*U(1,:).*U(2,:);
