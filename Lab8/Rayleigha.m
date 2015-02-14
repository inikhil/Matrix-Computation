function [iter,rho] = Rayleigh(A,q,k)

n = length(q);
iter = zeros(n,k);
sigma = zeros(k,1);
[Q,H] = hess(A);
iter(:,1) = q;
[m,r] = max(abs(q));
sigma(1) = q(r);
iter(:,1) = iter(:,1)/sigma(1);
rho = q'*(H*q);
rho = rho/(q'*q);
for j=1:(k-1)
    iter(:,j+1) = (H - (rho*eye(n)))\iter(:,j);
    [m,r] = max(abs(iter(:,j+1)));
    sigma(j+1) = iter(r,j+1);
    iter(:,j+1) = iter(:,j+1)/sigma(j+1);
    rho = iter(:,j+1)'*(H*iter(:,j+1))/(iter(:,j+1)'*iter(:,j+1));
end

[V,D] = eig(A);
d = diag(D);
[m,r] = max(abs(d));
v = V(:,r);
rate = norm(iter(:,k) - v,2)/norm(iter(:,k-1) - v,2);
disp(rate);

end