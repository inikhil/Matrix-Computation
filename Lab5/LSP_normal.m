function [x] = LSP_normal(t)

format long e;

m = length(t);
n = 20;

A = zeros(m,n);
count = 0;

for i=1:n
    A(:,i) = t.^count;
    count = count+1;
end

A1 = (A')*A;
b = sin(pi*t/5) + (t/5);
b1 = (A')*b';

disp(cond(A1,2));
x = A1\b1;

x = flipud(x);

err = polyval(x,t) - b;
disp(norm(err,2));

end