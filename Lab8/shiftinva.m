function [iter,sigma] = shiftinv(A,q,k,s)

n = length(q);
iter = zeros(n,k);
sigma = zeros(k,1);
[P,L,U] = lu(A-(s*eye(n)));
iter(:,1) = q;
[m,r] = max(abs(q));
sigma(1) = q(r);
iter(:,1) = iter(:,1)/sigma(1);
for j=1:(k-1)
    b = P*iter(:,j);
    y = forward_col_lower(L,b);
    iter(:,j+1) = backward_col_upper(U,y);
    [m,r] = max(abs(iter(:,j+1)));
    sigma(j+1) = iter(r,j+1);
    iter(:,j+1) = iter(:,j+1)/sigma(j+1);
end

sigma = (1./sigma) + s;
[V,D] = eig(A);
d = diag(D);
[m,r] = max(abs(d));
v = V(:,r);
rate = norm(iter(:,k) - v,2)/norm(iter(:,k-1) - v,2);
disp(rate);

end