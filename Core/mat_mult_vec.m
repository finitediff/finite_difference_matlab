function out = mat_mult_vec(A,B)

n = size(A,1);
N = size(B,2);

out = zeros(n,N);

for j = 1:n
    temp = zeros(1,N);
    for k = 1:n
        size(squeeze(A(j,k,:)))
        size(squeeze(B(k,:)))
        temp = temp+squeeze(A(j,k,:)).'.*squeeze(B(k,:));
    end
    out(j,:) = temp;
end