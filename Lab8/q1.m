A = [2.4347 -0.8641 1.4661;-0.8641 2.6630 -3.1263;1.4661 -3.1263 6.9023];
%A=[1 1 1; -1 9 2; -4 -1 2];
q = [1;1;1];
k = 200;
[iterate, sigma] = powermethod(A, q, k);
%norm(iterate(:,90),Inf)/norm(iterate(:,89),Inf)
%[V,D]=eig(A)
%abs(D(2,2))/abs(D(1,1))
[V,D] = eig(A);
d = diag(D);
[m,r] = max(abs(d));
v = V(:,r);
rate = norm(iterate(:,k) - v,Inf)/norm(iterate(:,k-1) - v,Inf);
disp(rate);
d1=sort(d);
d1(2)/d1(3)
