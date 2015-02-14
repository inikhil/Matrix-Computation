function [x] = LSP_aug(t)

format long e;

m = length(t);
n = 20;

A = zeros(m,n);
count = 0;

for i=1:n
    A(:,i) = t.^count;
    count = count+1;
end

A1(1:m,1:m) = eye(m);
A1(1:m,m+1:m+n) = A;
A1(m+1:m+n,1:m) = A';
A1(m+1:m+n,m+1:m+n) = zeros(n,n);
b = sin(pi*t/5) + (t/5);
b1(1:m,1) = b;
b1(m+1:m+n,1) = zeros(n,1);

x1 = A1\b1;

disp(cond(A1,2));

x = x1(m+1:m+n);

x = flipud(x);

err = polyval(x,t) - b;
disp(norm(err,2));

end