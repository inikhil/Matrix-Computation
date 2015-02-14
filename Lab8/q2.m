A = [2.4347 -0.8641 1.4661;-0.8641 2.6630 -3.1263;1.4661 -3.1263 6.9023];
q = [1;1;1];
k = 90;
shift = 100;
[iterate, sigma] = shiftinv(A, q, k, shift);
[V,D] = eig(A);
d = diag(D);
[m,r] = max(abs(d));
v = V(:,r);
rate = norm(iterate(:,k) - v,2)/norm(iterate(:,k-1) - v,2);
disp(rate);
d1=sort(d);
d1(2)/d1(3)